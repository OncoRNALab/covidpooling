### Load neccesary packages
library(tidyverse)

### Load and get data ready for the simulations
fw <- read_csv("../../Data/data_first_wave.csv")
sw <- read_csv("../../Data/data_second_wave.csv")

for_sim_fw <- fw %>%
  mutate(week = 1+abs(as.double(round(difftime(as.Date(strftime(min(first_wave$date), format="%Y-%m-%d")), date, units = "days") / 7)))) %>%
  group_by(week) %>%
  summarise(posrate = sum(cq != 0) / n(), 
            posval = list(cq[cq != 0]), 
            maxday = max(date), 
            quant = quantile(cq[cq != 0], p = 0.75)) %>%
  ungroup() %>%
  filter(week != 10)


for_sim_sw <- sw %>%
  mutate(week = 1+abs(as.double(round(difftime(as.Date(strftime(min(second_wave$date), format="%Y-%m-%d")), date, units = "days") / 7)))) %>%
  group_by(week) %>%
  summarise(posrate = sum(cq != 0) / n(), 
            posval = list(cq[cq != 0]), 
            maxday = max(date),
            quant = quantile(cq[cq != 0], p = 0.75)) %>%
  ungroup() %>%
  filter(week != 11)


### Define functions for simulations

#Simulate a dataset
create_data <- function(N, positive, cq_pos) {
  data <- tibble(index = 1:N)
  sample <- sample(data$index, positive)
  data <- data %>%
    mutate(posneg = ifelse(index %in% sample, 1, NA))
  data[!is.na(data$posneg), "posneg"] <- sample(cq_pos, sum(data$posneg, na.rm = TRUE), replace = TRUE)
  data
} 



#Create the groups (here indicated as plates)
create_plates <- function(data, plate_row, plate_col) {
  n <- nrow(data)
  s <- c()
  i <- 1
  while (n != 0) {
    if (n >= plate_row * plate_col) {
      n <- n - (plate_row * plate_col)
      s <- c(s, rep(i,plate_col * plate_row))
      i = i + 1
    } else {
      s <- c(s, rep(i, n))
      n <- 0
    }
  }
  
  data %>%
    mutate(plate = sample(s, length(s))) %>%
    group_by(plate) %>%
    mutate(rowcol = sample(1:n(), n())) %>%
    ungroup() %>%
    mutate(row = ceiling(rowcol / plate_col), 
           col = ifelse(rowcol %% plate_col == 0, plate_col, rowcol %% plate_col))
}

#Test sample with defined pooling approach
matrix_ntest <- function(plates, platesize, onrow = FALSE, oncol = FALSE) {
  ntest <- 0 #initialise the number of tests
  nplates <- 0 #initialise the number of plates used
  o_plates <- plates #remember orignal plates
  
  if (oncol) { #if ncol is true
    tested_plates_col <- plates %>% #group the columns on different plates and test them
      group_by(plate, col) %>%
      summarise(n = n(), l = sum(2^-posneg, na.rm = TRUE), post_cq = log2(n()) - log2(l),
                test = ifelse(post_cq >= 37 | is.na(post_cq), FALSE, TRUE)) %>%
      ungroup()
    ntest <- ntest + nrow(tested_plates_col) #add the number of tested columns to the amount of tests performed
    negative_col <- tested_plates_col %>% #find the columns on the plate that tested negative
      filter(!test) %>%
      dplyr::select(plate, col)
  }
  
  if (onrow) { #if onrow is true
    tested_plates_row <- plates %>%  #group the rows on different plates and test them
      group_by(plate, row) %>%
      summarise(n = n(), l = sum(2^-posneg, na.rm = TRUE), post_cq = log2(n()) - log2(l),
                test = ifelse(post_cq >= 37 | is.na(post_cq), FALSE, TRUE)) %>%
      ungroup() #test will be positive if any of the samples in the row was positive
    ntest <- ntest + nrow(tested_plates_row) #add the number of tested rows to the amount of tests performed
    negative_row <- tested_plates_row %>% #find the rows on the plate that tested negative
      filter(!test) %>%
      dplyr::select(plate, row)
  }
  
  if (oncol) {
    plates <- plates %>% anti_join(negative_col, by = c("plate" = "plate", "col" = "col"))
  } 
  if (onrow) {
    plates <- plates %>% anti_join(negative_row, by = c("plate" = "plate", "row" = "row"))
  }
  
  x <- o_plates %>%
    mutate(true = ifelse(posneg == 0, FALSE, TRUE)) %>%
    full_join(plates, by = c("index" = "index")) %>%
    mutate(tested = ifelse(is.na(posneg.y) | posneg.y == 0, FALSE, TRUE))
  falseneg <- 1 - mean(x %>%
                         filter(true == TRUE) %>%
                         .$tested)
  
  nplates <- ceiling(ntest / platesize)
  nplates <- nplates + ceiling(nrow(plates) / platesize)
  list(ntest = ntest + nrow(plates), #add the number of remaining samples to the number of test
       nplates = nplates,
       falseneg = falseneg,
       positive_values = x %>% filter(true == TRUE) %>% dplyr::select(posneg.x, tested, plate.x, row.x, col.x)) 
}

#Function to call for one simulation
sim <- function(N = 100000, pos = 5000, plate_col, plate_row, cq_pos) {
  s <- sapply(pos, function(x) {
    oncol <- ifelse(plate_row == 1, FALSE, TRUE)
    onrow <- ifelse(plate_col == 1, FALSE, TRUE)
    data <- create_data(N = N, positive = x, cq_pos)
    plates <- create_plates(data, plate_row, plate_col) 
    m <- matrix_ntest(plates, plate_row * plate_col, onrow = onrow, oncol = oncol)
    rbind(m, x) 
  })
  as_tibble(t(s)) %>%
    mutate(ntest = unlist(V1), nplate = unlist(V3), falseneg = unlist(V5)) %>%
    mutate(row = plate_row, col = plate_col, eff =  100000 / ntest) %>%
    dplyr::select(row, col, eff, falseneg)
}

#Function that does the simulation
simulations <- function(for_sim) {
  s <- tibble()
  for (i in 1:nrow(for_sim)) {
    set.seed(1)
    s1 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
    for (pool in pools[-1]) {
      set.seed(1)
      s_new <- sim(plate_row = pool[1], plate_col = pool[2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
      s1 <- rbind(s1, s_new)
    }
    
    set.seed(2)
    s2 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
    for (pool in pools[-1]) {
      set.seed(2)
      s_new <- sim(plate_row = pool[1], plate_col = pool[2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
      s2 <- rbind(s2, s_new)
    }
    
    set.seed(3)
    s3 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
    for (pool in pools[-1]) {
      set.seed(3)
      s_new <- sim(plate_row = pool[1], plate_col = pool[2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
      s3 <- rbind(s3, s_new)
    }
    
    set.seed(4)
    s4 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
    for (pool in pools[-1]) {
      set.seed(4)
      s_new <- sim(plate_row = pool[1], plate_col = pool[2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
      s4 <- rbind(s4, s_new)
    }
    
    set.seed(5)
    s5 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
    for (pool in pools[-1]) {
      set.seed(5)
      s_new <- sim(plate_row = pool[1], plate_col = pool[2], pos = 100000 * for_sim$posrate[i], cq_pos = unlist(for_sim$posval[i]))
      s5 <- rbind(s5, s_new)
    }
    
    #Combine everything
    s_new <- rbind(s1 %>% mutate(rep = 1), 
                   s2 %>% mutate(rep = 2),
                   s3 %>% mutate(rep = 3), 
                   s4 %>% mutate(rep = 4),
                   s5 %>% mutate(rep = 5)) %>%
      mutate(sensitivity = 1 - falseneg) %>%
      group_by(row, col) %>%
      summarise(eff = median(eff),
                min_eff = min(eff),
                max_eff = max(eff),
                sensitivity = median(sensitivity),
                min_sens = min(sensitivity),
                max_sens = max(sensitivity)) %>%
      ungroup() %>%
      mutate(week = for_sim$week[i])
    
    s <- rbind(s, s_new)
  }
  s
}



### Simulations
pools <- list(c(1,4), c(1,8), c(1, 12), c(1, 16), c(1, 24), c(8, 12), c(12, 16), c(16, 24)) #define pooling strategies

sim_fw <- simulations(for_sim = for_sim_fw) #simulate for the first dataset
sim_sw <- simulations(for_sim = for_sim_sw) #simulate for the second dataset






