library(haven)
library(dplyr)
library(forecast)
library(nlme)
library(psych)
library(ggplot2)
library(GGally)
library(leaps)
library(olsrr)
library(glmnet)
library(EnvStats)
#library(MASS)

data1 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data1_Srijan.dta")
data2 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data2_Srijan.dta")
data3 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data3_Srijan.dta")
data4 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data4_Srijan.dta")
data5 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data5_Srijan.dta")
data6 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data6_Srijan.dta")
data7 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data7_Srijan.dta")
data8 = read_dta("F:/B3_resources/ass_or_hw/Eco Stat/Data/Only_SrChattishgarh/data8_Srijan.dta")
data1_sr1 = na.omit(data1[, c(6, 19, 23, 28, 30, 31)])
data2_sr1 = na.omit(data2[, c(6, 19, 21, 25, 27, 28)])
data3_sr1 = na.omit(data3[, c(6, 20, 28, 32, 34, 35)])
data4_sr1 = na.omit(data4[, c(6, 19, 20, 24, 26, 27)])
data5_sr1 = na.omit(data5[, c(6, 19, 20, 24, 26, 27)])
summary = na.omit(data6[, c(6, 19, 20, 24, 26, 27)])

colname = colnames(data1_sr1)
colname[3] = "TCV"
colnames(data1_sr1) = colnames(data2_sr1) = colnames(data3_sr1) = colnames(data4_sr1) = colnames(data5_sr1) = colnames(summary) = colname


total = subset(summary, Item_Code == "40")

final_0 = bind_rows(data1_sr1, data2_sr1, data3_sr1, data4_sr1, data5_sr1)

repeatation = c(
  129,
  139,
  159,
  169,
  179,
  189,
  199,
  219,
  239,
  249,
  269,
  279,
  289,
  299,
  309,
  319,
  329,
  349,
  379,
  389,
  399,
  409,
  419,
  429,
  439,
  449,
  459,
  479,
  499,
  519,
  529,
  549,
  559,
  569,
  579,
  599,
  609,
  619,
  629,
  639,
  649,
  659
)

colnames(total)[3] <- "Total"


condition <- as.numeric(final_0$Item_Code) %in% repeatation

final = final_0[condition, ]

final <- final %>%
  left_join(dplyr::select(total, HHID, Total), by = "HHID")

final = final[, -1]

colnames(final)[2] = "TCV"
colnames(final)[6] = "Total"


engel_function = function(itemno) {
  item = subset(final, Item_Code == itemno)
  hh11 = subset(data7, HHID %in% sort(unique(as.numeric(item$HHID))))
  hh12 = subset(data8, HHID %in% sort(unique(as.numeric(item$HHID))))
  item <- item %>% arrange(HHID)
  hh11 <- hh11 %>% arrange(HHID)
  hh12 <- hh12 %>% arrange(HHID)
  item$size = hh11$hh_size
  item$earn = hh12$Regular_salary_earner
  item$mpce = hh12$MPCE
  item$mlt = hh12$MLT
  item$land = hh11$Land_Owned
  item$land[is.na(item$land)] <- 0
  df = cbind(
    as.numeric(item$TCV),
    as.numeric(item$Total),
    as.numeric(item$size),
    as.numeric(item$mpce),
    as.numeric(item$mlt),
    item$land
  )
  df = as.data.frame(df)
  colnames(df) <- c("TCV", "Total", "hhsize", "mpce", "mlt", "land")
  boxcox_result = EnvStats::boxcox(lm(TCV + 0.005 ~ ., data = df), optimize = T)
  lambda = boxcox_result$lambda
  df$TCV = boxcoxTransform(df$TCV + 0.005, lambda)
  degree1 = seq(0, 2, length = 100)
  degree2 = seq(0, 1, length = 5)
  degree3 = seq(0, 1, length = 3)
  r_squared_values <-
    array(0, dim = c(length(degree1), length(degree2), length(degree3)))
  for (i in 1:length(degree1)) {
    for (j in 1:length(degree2)) {
      for (k in 1:length(degree3)) {
        X_i <- df$Total ^ (degree1[i])
        Y_i <- df$hhsize ^ (degree2[j])
        Z_i <- df$land ^ (degree3[k])
        
        design_matrix <- cbind(1, X_i, Y_i, Z_i, df$mpce, df$mlt)
        
        model <- lm(df$TCV ~ ., data = as.data.frame(design_matrix))
        r_squared_values[i, j, k] <- summary(model)$r.squared
      }
    }
  }
  best_combi = which(r_squared_values == max(r_squared_values, na.rm = T), arr.ind = T)
  best = c(degree1[best_combi[1]], degree2[best_combi[2]], degree3[best_combi[3]])
  
  model1 <-
    lm(TCV ~ I((Total) ^ (degree1[best_combi[1]])) + I((hhsize) ^ (degree2[best_combi[2]])) +
         I((land) ^ (degree3[best_combi[3]])) + mpce + mlt, df)
  
  return_list1 = list(pairplot = ggpairs(df), Error = "Error!! Maximum R squared is also very less, A good model is not possible to fit!!")
  
  return_list2 = list(
    pairplot = ggpairs(df),
    best_model = model1,
    best_powers = best,
    boxcoxlambda = lambda
  )
  
  if (summary(model1)$r.squared < 0.5)
    return(return_list1)
  else
    (return(return_list2))
}




##Cereals

cereals = engel_function(129)
cereals$pairplot
x_limits <- c(min(df$Total), max(df$Total))
par(bg = 'lightgray')

densitytot = density(df$Total)
density_function = function(z) {
  densitytot$y[which.min(abs(densitytot$x - z))]
}

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 2),
  xlab = "Total Expenditure",
  ylab = "Expenditure Share for Cereals",
  main = "Engel Curves for Cereals"
)
f = function(x) {
  ((
    -169.940 + 88.442 * x ^ (0.1010101) + 2.759 * 2 + 0.067 * max(df$land) ^
      (0.5)
  ) ^ (1 / 0.5118896)) / x
}


curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)

for (i in hhseq) {
  
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((
      -169.940 + 88.442 * x ^ (0.1010101) + 2.759 * i + 0.067 * mean(df$land) ^
        (0.5)
    ) ^ (1 / 0.5118896)) / x
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
  
}


x_limits <- c(min(df$Total), max(df$Total))

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 20000),
  xlab = "Total Expenditure",
  ylab = "Expenditure for Cereals",
  main = "Engel Curves for Cereals"
)
f = function(x) {
  ((
    -169.940 + 88.442 * x ^ (0.1010101) + 2.759 * 2 + 0.067 * max(df$land) ^
      (0.5)
  ) ^ (1 / 0.5118896))
}
hhseq = seq(min(df$hhsize), max(df$hhsize), 1)
index_cereal = numeric(length(hhseq))
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)
k = 0
for (i in seq(min(df$hhsize), max(df$hhsize), 1)) {
  k = k + 1
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((
      -169.940 + 88.442 * x ^ (0.1010101) + 2.759 * i + 0.067 * mean(df$land) ^
        (0.5)
    ) ^ (1 / 0.5118896))
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
  seq <- seq(min(df$Total), max(df$Total), length = 1000)
  C1 = max(f(seq))
  g = function(x) {
    ((C1 - f(x)) / C1) * density_cerealfunction(x)
  }
  index_cereal[k] = sum(unlist(lapply(seq, g)))
}


##intoxicants

intoxicants = engel_function(329)



##Clothing

clothing = engel_function(379)

x_limits <- c(min(df$Total), max(df$Total))
par(bg = 'lightgray')

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 1.1),
  xlab = "Total Expenditure",
  ylab = "Expenditure Share for Clothing",
  main = "Engel Curves for Clothing"
)
f = function(x) {
  ((-410 + 42 * x ^ (0.6262626) - 1597 * 2 ^ (0.25) - 7 * max(df$land) ^ (0.5)) ^
     (1 / 1.000647)) / x
}
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)


for (i in seq(min(df$hhsize), max(df$hhsize), 1)) {
  
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((
      -410 + 42 * x ^ (0.6262626) - 1597 * i ^ (0.25) - 7 * max(df$land) ^ (0.5)
    ) ^ (1 / 1.000647)) / x
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
  
}


x_limits <- c(min(df$Total), max(df$Total))

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 60000),
  xlab = "Total Expenditure",
  ylab = "Expenditure for Clothing",
  main = "Engel Curves for Clothing"
)
f = function(x) {
  ((-410 + 42 * x ^ (0.6262626) - 1597 * 2 ^ (0.25) - 7 * max(df$land) ^ (0.5)) ^
     (1 / 1.000647))
}
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)
index_clothing = numeric(length(hhseq))
k = 0
for (i in hhseq) {
  k = k+1
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((
      -410 + 42 * x ^ (0.6262626) - 1597 * i ^ (0.25) - 7 * max(df$land) ^ (0.5)
    ) ^ (1 / 1.000647))
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
  seq <- seq(5000, max(df$Total), length = 1000)
  C1 = max(f(seq))
  g = function(x) {
    ((C1 - f(x)) / C1) * density_function(x)
  }
  index_clothing[k] = sum(unlist(lapply(seq, g)))
}



#clothing and cereals

index = as.data.frame(cbind(hhseq,index_cereal1,index_clothing1,index_dur1))
colnames(index)<-c("hhsize","cereals","clothing","durable")


##toilet

toilet = engel_function(459)


x_limits <- c(min(df$Total), max(df$Total))
par(bg = 'lightgray')

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 0.04),
  xlab = "Total Expenditure",
  ylab = "Expenditure Share for Toilet items",
  main = "Engel Curves for Toilet items"
)
f = function(x) {
  ((
    42.807 + 1.623 * x ^ (0.6060606) - 86.524 * 2 ^ (0.25) - 0.558 * max(df$land) ^
      (0.5)
  ) ^ (1 / 1.000647)) / x
}
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)
for (i in seq(min(df$hhsize), max(df$hhsize), 1)) {
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((
      42.807 + 1.623 * x ^ (0.6060606) - 86.524 * i ^ (0.25) - 0.558 * max(df$land) ^
        (0.5)
    ) ^ (1 / 1.000647)) / x
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
}


x_limits <- c(min(df$Total), max(df$Total))

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 1400),
  xlab = "Total Expenditure",
  ylab = "Expenditure for Toilet",
  main = "Engel Curves for Toilet"
)
f = function(x) {
  ((
    42.807 + 1.623 * x ^ (0.6060606) - 86.524 * 2 ^ (0.25) - 0.558 * max(df$land) ^
      (0.5)
  ) ^ (1 / 1.000647))
}
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)
index_dur = numeric(length(hhseq))
k = 0
for (i in seq(min(df$hhsize), max(df$hhsize), 1)) {
  k = k +1
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((
      42.807 + 1.623 * x ^ (0.6060606) - 86.524 * i ^ (0.25) - 0.558 * max(df$land) ^
        (0.5)
    ) ^ (1 / 1.000647))
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
  seq <- seq(5000, max(df$Total), length = 1000)
  C1 = max(f(seq))
  g = function(x) {
    ((C1 - f(x)) / C1) * density_function(x)
  }
  index_dur[k] = sum(unlist(lapply(seq, g)))
}



##Misc

misc = engel_function(499)

x_limits <- c(min(df$Total), max(df$Total))
par(bg = 'lightgray')

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 0.09),
  xlab = "Total Expenditure",
  ylab = "Expenditure Share for consumer services excluding conveyance items",
  main = "Engel Curves for consumer services excluding conveyance items"
)
f = function(x) {
  ((502 + x ^ (0.707) - 586 * 2 ^ (0.25)) ^ (1 / 1.000647)) / x
}
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)
for (i in seq(min(df$hhsize), max(df$hhsize), 1)) {
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((502 + x ^ (0.707) - 586 * i ^ (0.25)) ^ (1 / 1.000647)) / x
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
}


x_limits <- c(min(df$Total), max(df$Total))

plot(
  0,
  xlim = x_limits,
  ylim = c(0, 2600),
  xlab = "Total Expenditure",
  ylab = "Expenditure for consumer services excluding conveyance items",
  main = "Engel Curves for consumer services excluding conveyance items"
)
f = function(x) {
  ((502 + x ^ (0.707) - 586 * 2 ^ (0.25)) ^ (1 / 1.000647))
}
curve(f,
      from = min(df$Total),
      to = max(df$Total),
      add = T)
for (i in seq(min(df$hhsize), max(df$hhsize), 1)) {
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    ((502 + x ^ (0.707) - 586 * i ^ (0.25)) ^ (1 / 1.000647))
  }
  curve(
    f,
    min(df$Total),
    max(df$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
}




##durable

durable = engel_function(659)





data <- rnorm(100)
density_est <- density(data)
point_of_interest <- 0
density_at_point <-
  density_est$y[which.min(abs(density_est$x - point_of_interest))]
print(density_at_point)
















x_limits <- c(min(df1$Total), max(df1$Total))
y_limits <- c(30.3, 37)


plot(
  0,
  xlim = x_limits,
  ylim = y_limits,
  xlab = "X-axis",
  ylab = "Y-axis",
  main = "Engel Curves"
)


for (i in seq(min(df1$hhsize), max(df1$hhsize), 1)) {
  for (j in seq(min(df1$land), max(df1$land), length = 10)) {
    for (k in seq(min(df1$mlt), max(df1$mlt), length = 5)) {
      f = function(x) {
        x + x ^ (1 / 55) + i ^ (1 / 5) + j ^ (1 / 3) + k
      }
      curve(f, min(df1$Total), max(df1$Total), add = T)
    }
  }
}



x_limits <- c(min(df3$Total), max(df3$Total))
y_limits <- c(1.4, 1.8)
plot(
  0,
  xlim = x_limits,
  ylim = y_limits,
  xlab = "Total Expenditure",
  ylab = "Expenditure Share for Pulses",
  main = "Engel Curves for Pulses"
)

for (i in seq(min(df3$hhsize), max(df3$hhsize), 1)) {
  random_rgb <- sample(0:255, 3)
  random_color <-
    rgb(random_rgb[1], random_rgb[2], random_rgb[3], maxColorValue = 255)
  f = function(x) {
    0.518 * x + 0.006 * i + 0.436
  }
  curve(
    f,
    min(df3$Total),
    max(df3$Total),
    add = T,
    col = random_color,
    lwd = 2
  )
}