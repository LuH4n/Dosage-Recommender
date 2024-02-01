#Shiny dosage recommendation platform, written by Lu Han, 2024.01.10.
library(shiny)
library(rxode2)
library(foreach)
library(doParallel)
library(ggplot2)
library(lotri)
library(datasets)
library(dplyr)
library(tidyr)
library(data.table)

no_cores <- detectCores() - 1  # 使用的核心数，通常设置为总核心数减1
registerDoParallel(cores = no_cores)


#userinterface 用户交互界面界面设置
ui <- fluidPage(
  #标题设置
  titlePanel("Tacrolimus dosage recommendation system"),
  #侧边栏设置
  sidebarLayout(
    sidebarPanel(
      numericInput("weight", "Weight (kg)", value = 25),
      numericInput("rbc", "Red Blood Cell Count (10^12/mL)", value = 4.5),
      radioButtons("azole", "Azole Antifungal Agents", choices = c(" Co-administration" = 1, "Not Co-administration" = 0)),
      numericInput("therapeuticWindowLow", "Therapeutic Window Low (ng/mL)", value = 10),
      numericInput("therapeuticWindowHigh", "Therapeutic Window High (ng/mL)", value = 15),
      numericInput("dosingInterval", "Dosing Interval (h)", value = 12),
      actionButton("submit", "Submit")
    ),
    #主面板设置
    mainPanel(
      fluidRow(
        column(10, 
               wellPanel(
                 HTML('<h4>Recommended Dose:</h4>'),
                 textOutput("recommendedDose", container = span)
               )
        ),
        column(10, plotOutput("ctCurve"))
      )
    ) 
  )
)

#服务器设置
server <- function(input, output, session) {
  #自定义计算推荐剂量方程，需要用到从用户界面收集到的weight，rbc等信息。该方程的逻辑为尝试不同的治疗剂量，最终选定使得谷浓度最接近治疗窗的剂量
  calculateRecommendedDose <- function(weight, rbc, azole, therapeuticWindowLow, therapeuticWindowHigh, dosingInterval) {
    #定义合用啊唑类抗菌药物的影响，若合用则CL减52%
    azoleEffect <- ifelse(azole == "1", 0.522, 0) 
    #定义群体药代动力学参数，CL,V,KA
    POPpar <- c(14.3, 299, 4.48) 
    #使用rxode2包定义模型结构
    mod1 <- rxode2({
      CL = TVCL * (WT/26)^0.75 * (DDOSE/3)^0.358 * (RBC/3.26)^-1.71 * (1 - azoleEffect); 
      V2 = TVV2 * (WT/26)^1; 
      KA = KA; 
      S2 = V2/1000
      con = centr/V2
      conc = con*1000
      d/dt(depot1) = -KA*depot1; 
      d/dt(centr) = KA*depot1 - CL*con; 
    })
    #定义治疗靶标，靶标为治疗窗上下限之和除以二
    therapeuticWindowMid <- (therapeuticWindowLow + therapeuticWindowHigh) / 2
    #定义尝试的单次剂量，最低剂量为0.25，最高为12，步长为0.25
    doseOptions <- seq(0.25, 12, by = 0.25)
    #定义最接近的剂量这个变量，初始值为NULL
    closestDose <- NULL
    #定义尝试的谷浓度与治疗窗中点的差异这个变量，初始值为NULL
    minDifference <- Inf
    #循环尝试之前定义的单次剂量，得到谷浓度
    for (dose in doseOptions) {
      DD <- dose * 24 / dosingInterval #计算一天内的总剂量
      theta <- c(TVCL=POPpar[1], WT=weight, DDOSE=DD, RBC=rbc, azoleEffect=azoleEffect, 
                 TVV2=POPpar[2], 
                 KA=POPpar[3])
      #如果日剂量大于12则退出循环（一般对于儿童患者来说不建议使用大于12mg/d的他克莫司）
      if (DD > 12) {
        break
      }
      # 假设120小时的模拟时间（一般希望患者120达到治疗窗）
      dose_times <- seq(0, 120, by = dosingInterval) 
      # 设定谷浓度采样的时间（最后一次给药前0.001小时）
      specificTime <- 120 - 120 %% dosingInterval - 0.001
     # print(specificTime)print语句是用来debug的，如果结果有问题，可以输出中间过程的值，以发现问题
     #定义给药事件表格event table
      ev <- et(time = dose_times, amt = dose, cmt = 1)
      #在给药事件表格中加上采样事件，specificTime为预先定义的
      ev$add.sampling(specificTime)
      #print(ev)
      #使用solve函数进行模拟，模型为mod1，参数为theta，给药事件表格为ev，nSub=1代表仅模拟一次，由于这里没有应用个体间变异，因此使用均值进行模拟，只需要模拟一次
      sim1 <- solve(mod1, theta, ev, nSub=1)
      #print(sim1)
      #求得的谷浓度就是采样时间为specificTime的浓度
      troughConc <- sim1$conc[sim1$time == specificTime]
     
      #在这里，根据模型和输入参数计算 CL 和 V2
      calculatedCL <- POPpar[1] * (weight/26)^0.75 * (DD/3)^0.358 * (rbc/3.26)^-1.71 * (1 - azoleEffect)
      calculatedV2 <- POPpar[2] * (weight/26)
      
      #在这里计算刚刚提到的difference（各给药方案得到的谷浓度与治疗窗之间的差异）
      difference <- abs(troughConc - therapeuticWindowMid)
      #如果这一次循环得到的difference小于历次循环得到的minDifference，则minDifference替换为difference值，这次循环的dose就是最佳剂量方案
      if (difference < minDifference) {
        minDifference <- difference
        closestDose <- dose
      }
    }
    #函数的返回值即为bestDose等
    return(list(bestDose = closestDose, CL = calculatedCL, V2 = calculatedV2, AZOLE=azoleEffect))
  }
  
  #现在有了最佳的剂量bestDose（在这个函数里名称为optimalDose，但意思是一样的），使用bestDose做1000次模拟，以模拟患者可能的血药浓度范围，并生成药时曲线图
  simulateConcentration <- function(weight, rbc, azole, dosingInterval, optimalDose, nSim = 1000) {
    #这些代码的功能和上面calculateRecommendedDose函数代码的功能是一样的，只不过加了个体间变异eta1和eta2
    azoleEffect <- ifelse(azole == "1", 0.522, 0)
    POPpar <- c(14.3, 299, 4.48)
    CV_CL <- 0.387
    CV_V2 <- 0.875 
    DD <- optimalDose * 24 / dosingInterval
    
    mod <- rxode2({
      CL = TVCL * (WT/26)^0.75 * (DDOSE/3)^0.358 * (RBC/3.26)^-1.71 * (1 - zole)*exp(eta1); 
      V2 = TVV2 * (WT/26)^1 * exp(eta2); 
      KA = KA; 
      S2 = V2/1000
      con = centr/V2
      conc = con*exp(exp.err.pk)*1000
      d/dt(depot1) = -KA*depot1; 
      d/dt(centr) = KA*depot1 - CL*con; 
    })
    
    
    theta <- c(TVCL=POPpar[1], WT=weight, DDOSE=DD, RBC=rbc,zole =azoleEffect, TVV2=POPpar[2], KA=POPpar[3]) 
    #在这里设置个体间变异eta1和eta2值
    omega <- lotri(eta1~0.387^2,eta2~0.875^2)
    #模拟的残差设置为0，也可以在这里设置残差，
    sigma=lotri(exp.err.pk ~ 0)
    
    dose_times <- seq(0, 120, by = dosingInterval) # 假设120小时的模拟时间
    ev <- et(time = dose_times, amt = optimalDose, cmt = 1)
    #设置模拟的种子数
    set.seed(12345678)
    sim  <- solve(mod,theta, ev, omega=omega,sigma=sigma, nSub=1000,cores=no_cores)
    
    sim_data<- data.frame(sim)
    
    #收集用于绘制血药浓度曲线的数据，包括各个时间点的血药浓度点中位数，5%分位数和95%分位数
   result_plot <- sim_data %>%
      group_by(time) %>%
      summarise(
        c_0.05 = quantile(conc, probs = 0.05),
        c_0.5 = quantile(conc, probs = 0.5),
        c_0.95 = quantile(conc, probs = 0.95)
      )
    #这个函数返回的就是用于绘制血药浓度曲线的数据框
    return(result_plot)
  }
  
  #这个函数用于绘制血药浓度曲线
  plot_result <- function(result_plot,dosingInterval,therapeuticWindowLow,therapeuticWindowHigh) {
     #设置绘图的x轴最大值
    xMax <- (floor(120 / dosingInterval) * dosingInterval) + dosingInterval
     #用ggplot绘图
    ggplot(result_plot, aes(x = time)) + 
      geom_ribbon(aes(ymin = c_0.05, ymax = c_0.95), fill = "lightblue", alpha = 0.5) +  # 5%到95%的置信区间
      geom_line(aes(y = c_0.5), color = "blue") +  # 中位数曲线
      labs(x = "Time (days)", y = "Tacrolimus concentration (ng/mL)") + 
      scale_x_continuous(limits = c(0, xMax),breaks=c(0,24,48,72,96,120),labels = c(0,1,2,3,4,5)) +
      scale_y_continuous(breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70)) +
      geom_hline(yintercept = c(therapeuticWindowLow,therapeuticWindowHigh), linetype = 2, size = 0.6, colour = "skyblue")+ 
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", size = 12, color = "black", hjust = 0),
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text = element_text(face = "bold", size = 10, color = "black"),
        strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "grey"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major.y = element_line(color = "white", linetype = 2),
        panel.grid.major.x = element_line(color = "white", linetype = 2),
        panel.grid.minor.y = element_line(color = "white", linetype = 2),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.6, 0.6), "lines")
      )
  }
  
  #监听函数，当用户点击submit键时，从用户交互界面中获取值，并输出值
  observeEvent(input$submit, {
    #当用户点击submit键时，调用recommendedDose函数，得到的最优剂量存储在recommendedDose中
    recommendedDose <- calculateRecommendedDose(
      input$weight, 
      input$rbc, 
      input$azole, 
      input$therapeuticWindowLow, 
      input$therapeuticWindowHigh, 
      input$dosingInterval
    )
    #当用户点击submit键时，在用户交互界面输出recommendedDose，renderText函数用于渲染输出的文本
    output$recommendedDose <- renderText({
      paste0(recommendedDose$bestDose, " mg\n"," q",input$dosingInterval,"h")
    })
    #当用户点击submit键时，调用simulateConcentration函数和plot_result函数，在用户交互界面输出画出的药时曲线，renderPlot函数用于渲染输出的图片
    output$ctCurve <- renderPlot({
      result_plot <- simulateConcentration(input$weight, input$rbc, input$azole, input$dosingInterval, recommendedDose$bestDose)
      plot_result(result_plot,input$dosingInterval,input$therapeuticWindowLow, input$therapeuticWindowHigh)
    })
})
}
  
    # 运行Shiny应用程序
    shinyApp(ui = ui, server = server)
    
    #closeAllConnections()，如果报错All connection is used, 运行这个语句，再尝试重启R
    