#Load all necessary packages
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(rlist)
library(data.table)
library(zoo)
library(ggbeeswarm)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
gradient <- c("#5D5D8D", "#A2699E", "#E07897", "#FF9580", "#FFC36A", "#F9F871")
gradient <- c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080")

#Load data
data <- read_table2("../../Data/data_COVID_pooling2.txt", col_names = FALSE)
wells <- read_table2("../../Data/well_COVID_pooling2.txt", col_names = FALSE)
plates <- read_table2("../../Data/listplates.txt", col_names = FALSE)
counts <- read_table2("../../Data/counter.txt", col_names = FALSE)

s <- c()
for (i in 1:nrow(counts)) {
  s <- c(s, rep(plates$X1[i], counts$X1[i]))
}
x <- rep(s, each = 8)
x2 <- rep(1:length(s), each = 8)

new_data <- data %>%
  mutate(PCRplate = x, RNAplate = x2)
new_wells <- wells %>%
  mutate(PCRplate = x, RNAplate = x2)



d <- new_data %>%
  mutate(row = rep(1:8, nrow(new_data) / 8), #add a row number
         group = rep(1:(nrow(new_data) / 16), each = 16)) %>% #add to which PCR-plate the rows belong
  gather(key = "column", value = "cq", X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12) %>% #gather all columns 
  mutate(id = 1:n(), cq = ifelse(cq == 0, NA, cq)) #id all wells

w <- new_wells %>%
  mutate(row = rep(1:8, nrow(data) / 8), #add to which PCR-plate the rows belong
         group = rep(1:(nrow(data) / 16), each = 16)) %>%  #add a row number
  gather(key = "column", value = "wells", X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12) %>% #gather all columns 
  mutate(id = 1:n()) #id all wells

all_data <- d %>%
  full_join(dplyr::select(w, id, wells), by = c("id" = "id"))

nonfull_plate <- unique(all_data %>% #get all the plates that were not full
                          filter(is.na(wells)) %>%
                          .$RNAplate)

bad_control <- unique(all_data %>% #get all plates that have a bad control value
  filter((row == 8 & column == "X12" & is.na(cq)) | (row == 7 & column == "X12" & !is.na(cq))) %>%
  .$PCRplate)

con_ids <- all_data %>% #get the ids that are control wells
  filter(row == 8 & column == "X12" | row == 7 & column == "X12") %>%
  .$id 

many_plate <- all_data %>% #get plates with too many positives
  group_by(RNAplate) %>%
  summarise(pos = sum(!is.na(cq))) %>%
  ungroup() %>%
  filter(pos >= 10) %>%
  .$RNAplate

avg_all <- mean(all_data %>%
                  filter(!(RNAplate %in% many_plate), !(RNAplate %in% nonfull_plate),  !(RNAplate %in% bad_control), row == 8, column == "X12") %>%
                  .$cq, na.rm = TRUE)

correction <- all_data %>% 
  filter(row == 8 & column == "X12") %>%
  group_by(PCRplate) %>% #calculate the mean of the positive controls on each PCR-plate
  summarize(avg = mean(cq, na.rm = TRUE),diff = avg - avg_all) %>% #find how much the values in this plate should be adjusted
  ungroup()


min = mean(correction$avg, na.rm = TRUE) - 2*sd(correction$avg, na.rm = TRUE)
max = mean(correction$avg, na.rm = TRUE) + 2*sd(correction$avg, na.rm = TRUE)

good_cont <- correction %>% 
  filter(avg > min, avg < max) %>%
  .$PCRplate
#Filter out all plates with have too many positives, are not full or have bad control values. Also remove any control wels
d_fil <- all_data %>%
  filter(!(RNAplate %in% many_plate), !(id %in% con_ids), !(RNAplate %in% nonfull_plate),  !(RNAplate %in% bad_control), PCRplate %in% good_cont) %>%
  left_join(correction, by = c("PCRplate" = "PCRplate")) %>% #join the correction
  mutate(corr_cq = cq - diff) %>% #calculate the correction
  arrange(PCRplate, RNAplate, row) #just to make it nicer

d_fil_hosp <- all_data %>%
  filter(!(id %in% con_ids), !(RNAplate %in% nonfull_plate),  !(RNAplate %in% bad_control), PCRplate %in% good_cont) %>%
  left_join(correction, by = c("PCRplate" = "PCRplate")) %>% #join the correction
  mutate(corr_cq = cq - diff) %>% #calculate the correction
  arrange(PCRplate, RNAplate, row) #just to make it nicer

d_fil_ohosp <- all_data %>%
  filter((RNAplate %in% many_plate), !(id %in% con_ids), !(RNAplate %in% nonfull_plate),  !(RNAplate %in% bad_control), PCRplate %in% good_cont) %>%
  left_join(correction, by = c("PCRplate" = "PCRplate")) %>% #join the correction
  mutate(corr_cq = cq - diff) %>% #calculate the correction
  arrange(PCRplate, RNAplate, row) #just to make it nicer

#Get the corrected cq values. Since 37 is the cut-off, we will only use values under 37
cq_pos <- d_fil %>%
  filter(corr_cq <= 37) %>%
  .$corr_cq


data_pub <- d_fil %>%
  filter(corr_cq <= 37)  %>%
  mutate(target = "E gene", sample_type = "std", plate_id = PCRplate,
         target_type = "toi", sample = id, dye = "FAM", reaction = ifelse(plate_id %in% plates$X1[1:155], "singleplex", "duplex")) %>%
  dplyr::select(plate_id, sample, sample_type, target, target_type, dye, reaction, wells, corr_cq)
  

write_csv(data_pub, "../../Data/filtered_cq.csv")


d_cool <- rbind(d_fil, d_fil_hosp, d_fil_ohosp) %>%
  mutate(group = c(rep("without high prevalence", nrow(d_fil)), rep("with high prevalence", nrow(d_fil_hosp)), rep("only high prevalence", nrow(d_fil_ohosp))))
### plot for Jim

d_cool %>%
  ggplot(aes(corr_cq, col = group)) +
  stat_ecdf() +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 15)) +
  ylab("density\n") +
  xlab("\nCq value") 
  


d_fil %>%
  ggplot(aes(corr_cq)) +
  geom_histogram(fill = "#2D5580", binwidth = 1) +
  xlab("\nCq value") +
  ylab("frequency\n") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +
  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  ggtitle("Distribution for filtered Cq values", 
          subtitle = "filtered for only full 96-well RNA plates with good controls \nand less than 10 positives (N = 1676)") 


d_fil_hosp %>%
  ggplot(aes(corr_cq)) +
  geom_histogram(fill = "#2D5580", binwidth = 1) +
  xlab("\nCq value") +
  ylab("frequency\n") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +
  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  ggtitle("Distribution for filtered Cq values", 
          subtitle = "filtered for only full 96-well RNA plates with good controls (N = 4103)") 

####### Define important functions ##########
#Simulate a dataset
create_data <- function(N, positive) {
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
data = create_data(100, 10)
plates <- create_plates(data, 8, 12)
###### TEST THE SAMPLES MATRIX-WISE, EITHER BY COLUMN OR ROW, OR BOTH #######
matrix_ntest <- function(plates, platesize, onrow = FALSE, oncol = FALSE) {
  ntest <- 0 #initialise the number of tests
  nplates <- 0 #initialise the number of plates used
  o_plates <- plates #remember orignal plates
  
  if (oncol) { #if ncol is true
    tested_plates_col <- plates %>% #group the columns on different plates and test them
      mutate(egene = ifelse(is.na(posneg), 0, 2 ^ (40 - posneg))) %>% #group the rows on different plates and test them
      group_by(plate, col) %>%
      summarise(post_cq = 40 - log2(mean(egene)),
                test = ifelse(post_cq >= 37 | is.na(post_cq), FALSE, TRUE)) %>%
      ungroup()
    ntest <- ntest + nrow(tested_plates_col) #add the number of tested columns to the amount of tests performed
    negative_col <- tested_plates_col %>% #find the columns on the plate that tested negative
      filter(!test) %>%
      select(plate, col)
  }
  
  if (onrow) { #if onrow is true
    tested_plates_row <- plates %>% 
      mutate(egene = ifelse(is.na(posneg), 0, 2 ^ (40 - posneg))) %>% #group the rows on different plates and test them
      group_by(plate, row) %>%
      summarise(post_cq = 40 - log2(mean(egene)),
                test = ifelse(post_cq >= 37 | is.na(post_cq), FALSE, TRUE)) %>%
      ungroup() #test will be positive if any of the samples in the row was positive
    ntest <- ntest + nrow(tested_plates_row) #add the number of tested rows to the amount of tests performed
    negative_row <- tested_plates_row %>% #find the rows on the plate that tested negative
      filter(!test) %>%
      select(plate, row)
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
      select(plate, col)
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
      select(plate, row)
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
############## Simulation ###########

sim <- function(N = 100000, pos = round(10 ^ seq(-4, -1, by = 0.05) * 100000), plate_col, plate_row) {
  s <- sapply(pos, function(x) {
    oncol <- ifelse(plate_row == 1, FALSE, TRUE)
    onrow <- ifelse(plate_col == 1, FALSE, TRUE)
    data <- create_data(N = N, positive = x)
    plates <- create_plates(data, plate_row, plate_col) 
    m <- matrix_ntest(plates, plate_row * plate_col, onrow = onrow, oncol = oncol)
    rbind(m, x) 
  })
  as_tibble(t(s)) %>%
    mutate(ntest = unlist(V1), nplate = unlist(V3), falseneg = unlist(V5), pos = unlist(V4), positive_values = V7) %>%
    mutate(row = plate_row, col = plate_col) %>%
    dplyr::select(ntest, nplate, falseneg, pos, positive_values, row, col)
}  


pools <- list(c(1,4), c(1,8), c(1, 12), c(1, 16), c(1, 24), c(8, 12), c(12, 16), c(16, 24))

set.seed(1)
s1 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2])
for (pool in pools[-1]) {
  set.seed(1)
  s_new <- sim(plate_row = pool[1], plate_col = pool[2])
  s1 <- rbind(s1, s_new)
}

set.seed(2)
s2 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2])
for (pool in pools[-1]) {
  set.seed(2)
  s_new <- sim(plate_row = pool[1], plate_col = pool[2])
  s2 <- rbind(s2, s_new)
}

set.seed(3)
s3 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2])
for (pool in pools[-1]) {
  set.seed(3)
  s_new <- sim(plate_row = pool[1], plate_col = pool[2])
  s3 <- rbind(s3, s_new)
}

set.seed(4)
s4 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2])
for (pool in pools[-1]) {
  set.seed(4)
  s_new <- sim(plate_row = pool[1], plate_col = pool[2])
  s4 <- rbind(s4, s_new)
}

set.seed(5)
s5 <- sim(plate_row = pools[[1]][1], plate_col = pools[[1]][2])
for (pool in pools[-1]) {
  set.seed(5)
  s_new <- sim(plate_row = pool[1], plate_col = pool[2])
  s5 <- rbind(s5, s_new)
}


s <- rbind(s1 %>% mutate(rep = 1), 
           s2 %>% mutate(rep = 2),
           s3 %>% mutate(rep = 3), 
           s4 %>% mutate(rep = 4),
           s5 %>% mutate(rep = 5))


s1 %>%
  mutate(strategy = paste0(row, "x", col)) %>%
  group_by(strategy) %>%
  summarise(sum(pos))


positive <- rbind(list.rbind(s1$positive_values),
      list.rbind(s2$positive_values),
      list.rbind(s3$positive_values),
      list.rbind(s4$positive_values),
      list.rbind(s5$positive_values)) 


test <- c()
for (i in 1:nrow(s)) {
  test <- c(test, rep(s$pos[i], s$pos[i]))
}

strategies <- c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")
p <- positive %>%
  mutate(rep = rep(c("1", "2", "3", "4", "5"), each = nrow(positive) / 5),
         strategy = factor(rep(rep(strategies, each = 91871), 5), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         prevalence = test / 1000) 




s <- s %>%
  mutate(sensitivity = 1 - falseneg)



range_cq <- p %>%
  filter(tested == FALSE) %>%
  group_by(strategy) %>%
  summarise(min = min(posneg.x), max = max(posneg.x)) %>%
  ungroup() 



ja <- p %>%
  group_by(strategy, rep, posneg.x) %>%
  summarize(m = mean(tested)) %>%
  ungroup() %>%
  group_by(strategy, posneg.x) %>%
  summarise(avg = mean(m), min = min(m), max = max(m)) %>%
  mutate(rollmean = frollmean(avg, 20))

ja2 <- p %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(m = mean(tested)) %>%
  ungroup()

p1 <- p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  mutate(tested = as.factor(tested))
levels(p1$tested) <- c("FN", "TP")

p1 %>%
  ggplot(aes(posneg.x, linetype = tested, col = strategy)) +
  stat_ecdf() +
  facet_grid(strategy~prevalence) +
  theme_minimal() +
  ylab("cummultative distribution\n") +
  xlab("\nCq") +
  ggtitle("Distribution of Cqs > 32")



p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(TP = sum(tested == "TP"), total = n()) %>%
  ungroup() %>%
  group_by(strategy, prevalence) %>%
  mutate(cs = cumsum(TP), cs_total = cumsum(total), sens = cs / cs_total) %>%
  ggplot(aes(posneg.x, sens)) +
  geom_line() +
  facet_grid(strategy~prevalence) +
  theme_minimal() +
  ylab("cummultative sensitivity\n") +
  xlab("\nCq") +
  xlim(30, 37) +
  theme(axis.title = element_text(face = "bold", size  = 15),
        strip.text = element_text(face = "bold", size = 12))
  

p2 <- p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  group_by(strategy, prevalence) %>%
  summarize(n = n())


ggplot(aes(posneg.x, FN, col = strategy)) +
  geom_line() +
  facet_grid(strategy~prevalence) +
  theme_minimal() +
  ylab("cummultative distribution\n") +
  xlab("\nCq") +
  ggtitle("Distribution of Cqs > 32")

ja2 %>%
  filter(prevalence %in% c(0.01, 0.1, 1, 10)) %>%
  mutate(r = frollmean(m, 5)) %>%
  ggplot(aes(posneg.x, r, col = strategy)) +
  geom_line(linetype= "dashed") +
  xlim(30, 37) +
  facet_wrap(~prevalence, nrow = 2)
 
ja %>%
  ggplot(aes(posneg.x, rollmean)) +
  geom_line() +
  facet_grid(~strategy) +
  xlim(30, 37) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 15)) +
  ylab("rolling mean sensitivity\n") +
  xlab("\nCq")







t %>% 
  group_by(posneg.x) %>%
  summarise(m = mean(tested)) %>%
  ungroup() %>%
  mutate(avg = frollmean(m, 20)) %>%
  ggplot(aes(posneg.x, avg)) +
  geom_line(size = 1) +
  xlim(30, 37) +
  theme_minimal() +
  ylab("Rolling mean sensitivity (20 values)\n") +
  xlab("\nCq value") +
  theme(axis.title = element_text(face = "bold", size = 15))



a <- t %>% 
  group_by(posneg.x) %>%
  summarise(m = mean(tested)) %>%
  ungroup() %>%
  mutate(avg = frollmean(m, 50))


t %>%
  mutate(pool = posneg.x > 34) %>%
  filter(pool, tested) %>%
  mutate()

rowt <- t %>% group_by(prevalence, plate.x, row.x) %>%
  summarize(l = list(posneg.x)) %>%
  ungroup() 

q <- c()
for (i in 1:nrow(rowt)) {
  q <- c(q, any(unlist(rowt$l[i]) < 34))
}

rowt <- rowt %>%
  mutate(check = q)


t %>%
  inner_join(select(rowt, prevalence, plate.x, row.x, check),
             by = c("prevalence" = "prevalence", "row.x" = "row.x", "plate.x" = "plate.x")) %>%
  mutate(if_alone = posneg.x < 34) %>%
  filter(!if_alone) %>%
  group_by(tested, check) %>%
  summarize(n = n())



############ Figures ###############
s %>%
  mutate(eff = 100000 / ntest, group = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         Prevalence = pos / 100000 * 100,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  ggplot(aes(eff, sensitivity, col = Prevalence)) +
  geom_point() +
  facet_wrap(~group,
             nrow = 2) +
  ylab("sensitivity\n") +
  xlab("\nefficiency increase") +
  labs(col = "prevalence (%)\n") +
  theme_minimal() +
  scale_color_gradientn(trans = "log",
                        breaks = c(10, 1, 0.1, 0.01),
                        colours = gradient) +
  theme(legend.position = c(1.2, 0.2),
        legend.title = element_text(face = "bold",
                                    size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold",
                                  size = 15),
        axis.title = element_text(face = "bold",
                                  size = 17),
        axis.text = element_text(size = 15))
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure4A.png")




###### With middle sensitivity ############
zone_mean <- s %>%
  mutate(group = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         prevalence = pos / 100000 * 100) %>%
  filter(prevalence >= 0.01, prevalence <= 1) %>%
  group_by(group) %>%
  summarize(avg = median(sensitivity))

s %>%
  mutate(efficiency = 100000 / ntest, group = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         Prevalence = pos / 100000 * 100,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  group_by(group, Prevalence) %>%
  summarise(sens = median(sensitivity), eff = median(efficiency)) %>%
  ungroup() %>%
  ggplot(aes(Prevalence, sens, col = eff)) +
  geom_point() +
  geom_segment(data = zone_mean, aes(x = 0.1, xend = 1, y = avg, yend = avg), col = "black", size = 1) +
  geom_text(data = zone_mean, aes(x = 10^-0.5, y = avg + 0.05, label = round(avg, digits = 2)), col = "black", fontface = "bold") +
  facet_wrap(~group,
             nrow = 2,
             scales = "free_x") +
  ylab("sensitivity\n") +
  xlab("\nprevalence (%)") +
  labs(col = "efficiency") +
  theme_minimal() +
  scale_x_log10() +
  scale_color_gradientn(colours = gradient) +
  theme(legend.position = c(1.2, 0.5),
        legend.title = element_text(face = "bold",
                                    size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold",
                                  size = 15),
        axis.title = element_text(face = "bold",
                                  size = 17),
        axis.text = element_text(size = 8),
        panel.spacing = unit(1, "lines"),
        plot.margin=unit(c(0,4.5,0,1),"cm")) 

######## With linear regression ############@
s %>%
  mutate(efficiency = 100000 / ntest, group = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         Prevalence = pos / 100000 * 100,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  group_by(group, Prevalence) %>%
  summarise(sens = median(sensitivity), eff = median(efficiency)) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(intercept = lm(sens ~ Prevalence)$coefficients[1], 
         slope = lm(sens ~ Prevalence)$coefficients[2],
         pred = intercept + Prevalence * slope) %>%
  ungroup() %>%
  ggplot(aes(Prevalence, sens, col = eff)) +
  geom_point() +
  geom_line(aes(x = Prevalence, y = pred), col ="black") +
  geom_text(aes(label = paste0(round(intercept, digits = 2), " + ", round(slope,digits = 3), "p"), y = 0.55, x = 5),
            col = "black",
            size = 4) +
  facet_wrap(~group,
             nrow = 2,
             scales = "free_x") +
  ylab("sensitivity\n") +
  xlab("\nprevalence (%)") +
  ylim(0.5, 1.02) +
  labs(col = "efficiency") +
  theme_minimal() +
  scale_color_gradientn(colours = gradient) +
  theme(legend.position = c(1.2, 0.5),
        legend.title = element_text(face = "bold",
                                    size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold",
                                  size = 15),
        axis.title = element_text(face = "bold",
                                  size = 17),
        axis.text = element_text(size = 8),
        panel.spacing = unit(1, "lines"),
        plot.margin=unit(c(0,4.5,0,1),"cm")) 
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure4B2.png")

### Supplemental 3 #####
s %>%
  mutate(efficiency = 100000 / ntest, group = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         Prevalence = pos / 100000 * 100,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  group_by(group, Prevalence) %>%
  summarise(sens = median(sensitivity), eff = median(efficiency)) %>%
  ungroup() %>%
  ggplot(aes(group, sens, col = Prevalence)) +
  geom_beeswarm() +
  ylab("sensitivity\n") +
  xlab("\nprevalence (%)") +
  ylim(0.5, 1.02) +
  labs(col = "prevalence") +
  theme_minimal() +
  scale_color_gradientn(colours = gradient) +
  theme(legend.position = c(1.2, 0.5),
        legend.title = element_text(face = "bold",
                                    size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold",
                                  size = 15),
        axis.title = element_text(face = "bold",
                                  size = 17),
        axis.text = element_text(size = 8),
        panel.spacing = unit(1, "lines"),
        plot.margin=unit(c(0,4.5,0,1),"cm")) 


d_fil %>%
  ggplot(aes(corr_cq)) +
  geom_histogram(col = "#2D5580", binwidth = 1, fill = 'white', size = 1) +
  xlab("\nCq value") +
  ylab("frequency\n") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  scale_y_continuous(breaks = seq(0, 125, by = 25))
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure1.png")


s %>%
  mutate(strategy = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         prevalence = pos / 100000 * 100,
         efficiency = 100000 / ntest,
         plate_efficiency = 261 / nplate,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  group_by(dimension,strategy, prevalence) %>%
  summarise(avg = median(efficiency),
            sd = sd(efficiency),
            n = n(),
            ci = 1.96 * sd / sqrt(n),
            upper = avg + ci,
            lower = avg - ci,
            min = min(efficiency), 
            max = max(efficiency)) %>%
  ggplot(aes(prevalence, avg, col = strategy)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = strategy), alpha = 0.5, colour = NA) +
  geom_line(size = 0.5) +
  scale_x_log10() +
  ylab("average efficiency increase\n") +
  xlab("\nprevalence (%)") +
  theme_minimal() +
  theme(legend.position = "left",
        legend.title = element_text(face = "bold",
                                    size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.title = element_text(face = "bold",
                                  size = 17),
        axis.text = element_text(size = 15)) +
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  scale_linetype_manual(values=c(2,1)) 
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure2_2.png")

p <- p %>%
  mutate(tested = as.factor(tested)) 
levels(p$tested) = c("FN", "TP")


p %>%
  ggplot(aes(factor(as.character(tested), levels = c("FN", "TP")), posneg.x, col = strategy)) +
  geom_boxplot() +
  facet_grid(~strategy) +
  theme_minimal() +
  ylab("Cq\n") +
  theme(axis.title.y = element_text(face = "bold", size = 17),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values=cbPalette) 
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure5_2.png",
       width = unit(8, "cm"),
       height = unit(4, "cm"))



ja %>%
  ggplot(aes(posneg.x, rollmean, col = strategy)) +
  geom_line(size = 1) + 
  xlim(30, 37) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 17),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold")) +
  ylab("rolling mean sensitivity\n") +
  xlab("\nCq") +
  scale_colour_manual(values=cbPalette)
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure6.png")


p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(TP = sum(tested == "TP"), total = n()) %>%
  ungroup() %>%
  group_by(strategy, prevalence) %>%
  mutate(cs = cumsum(TP), cs_total = cumsum(total), sens = cs / cs_total) %>%
  ggplot(aes(posneg.x, sens)) +
  geom_line() +
  facet_grid(strategy~prevalence) +
  theme_minimal() +
  ylab("cummultative sensitivity\n") +
  xlab("\nCq") +
  xlim(30, 37) +
  theme(axis.title = element_text(face = "bold", size  = 15),
        strip.text = element_text(face = "bold", size = 10))
ggsave("/Users/jverwilt/Dropbox (OncoRNALab)/Basecamp/JAV/JAV2003 COVID Pooling/Report/figures/figure7_2.png")

######## With colours instead of faceted ########@
p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(TP = sum(tested == "TP"), total = n()) %>%
  ungroup() %>%
  group_by(strategy, prevalence) %>%
  mutate(cs = cumsum(TP), cs_total = cumsum(total), sens = cs / cs_total) %>%
  ggplot(aes(posneg.x, sens, col = strategy)) +
  geom_line() +
  facet_grid(~prevalence) +
  theme_minimal() +
  ylab("cummultative sensitivity\n") +
  xlab("\nCq") +
  xlim(30, 37) +
  theme(axis.title = element_text(face = "bold", size  = 15),
        strip.text = element_text(face = "bold", size = 10)) +
  scale_colour_manual(values=cbPalette)


p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  mutate(prevalence = factor(prevalence, levels = c("0.01", "0.1", "1", "10"))) %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(TP = sum(tested == "TP"), total = n()) %>%
  ungroup() %>%
  group_by(strategy, prevalence) %>%
  mutate(cs = cumsum(TP), cs_total = cumsum(total), sens = cs / cs_total) %>%
  ggplot(aes(posneg.x, sens, col = prevalence)) +
  geom_line() +
  facet_wrap(~strategy,
             nrow = 2) +
  theme_minimal() +
  ylab("cummultative sensitivity\n") +
  xlab("\nCq") +
  xlim(30, 37) +
  theme(axis.title = element_text(face = "bold", size  = 15),
        strip.text = element_text(face = "bold", size = 10)) +
  scale_colour_manual(values=cbPalette[])




hu <- ja2 %>%
  filter(prevalence %in% c(0.01, 0.1, 1, 10)) %>%
  group_by(strategy, prevalence) %>% 
  arrange(prevalence, strategy, posneg.x) %>%
  mutate(rollmean = frollmean(m, 5))


p %>%
  filter(prevalence %in% c(0.01, 0.1, 1, 10)) %>%
  arrange(strategy, prevalence, posneg.x) %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(total = n(), positives = sum(tested == "TP"), 
            ) %>%


hu %>%
  ungroup() %>%
  ggplot(aes(posneg.x, rollmean, col = strategy)) +
  geom_line() +
  facet_wrap(~prevalence, nrow = 2) +
  xlim(30, 37)



ja %>%
  ggplot(aes(posneg.x, rollmean)) +
  geom_line(size = 1) + 
  geom_point(data = ja2, aes(posneg.x, m)) +
  facet_wrap(~strategy, nrow = 2) +
  xlim(30, 37) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 17),
        axis.text = element_text(size = 15),
        legend.position = "none",
        strip.text = element_text(size = 15, face = "bold")) +
  ylab("rolling mean sensitivity\n") +
  xlab("\nCq") +
  scale_colour_manual(values=cbPalette)



what <- s %>%
  mutate(strategy = factor(paste0(row, "x", col), levels = c("1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         prevalence = pos / 100000 * 100,
         efficiency = 100000 / ntest,
         plate_efficiency = 261 / nplate,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  group_by(dimension,strategy, prevalence) %>%
  summarise(avg = median(efficiency),
            sd = sd(efficiency),
            n = n(),
            ci = 1.96 * sd / sqrt(n),
            upper = avg + ci,
            lower = avg - ci,
            min = min(efficiency), 
            max = max(efficiency))%>%
  ungroup() %>%
  group_by(prevalence) %>%
  summarize(maximum = max(avg), max_strat = strategy[which(avg == maximum)])

x <- c()
for (i in 1:2) {
  set.seed(i)
  x[i] = rnorm(1)
}

x <- c()
set.seed(1)
x[1] = rnorm(1)
set.seed(2)
x[2] = rnorm(1)
