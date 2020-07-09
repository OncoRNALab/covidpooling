#Load all necessary packages
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(rlist)
library(data.table)
library(zoo)
library(grid)
library(gtable)
library(ggbeeswarm)

#Define color palettes
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
gradient <- c("#5D5D8D", "#A2699E", "#E07897", "#FF9580", "#FFC36A", "#F9F871")
gradient <- c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080")

############## Generate starting dataset ##################@@
#Load data
data <- read_table2("../../Data/data_COVID_pooling2.txt", col_names = FALSE)
wells <- read_table2("../../Data/well_COVID_pooling2.txt", col_names = FALSE)
plates <- read_table2("../../Data/listplates.txt", col_names = FALSE)
counts <- read_table2("../../Data/counter.txt", col_names = FALSE)

#Data manipulation
s <- c()
for (i in 1:nrow(counts)) {
  s <- c(s, rep(plates$X1[i], counts$X1[i]))
}


d <- data %>%
  mutate(PCRplate = rep(s, each = 8), #add PCR plate information
         RNAplate = rep(1:length(s), each = 8), #add RNA plate information
         row = rep(1:8, nrow(data) / 8), #add a row number
         group = rep(1:(nrow(data) / 16), each = 16)) %>% #add to which PCR-plate the rows belong
  gather(key = "column", value = "cq", X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12) %>% #gather all columns 
  mutate(id = 1:n(), cq = ifelse(cq == 0, NA, cq)) #id all wells

w <- wells %>%
  mutate(PCRplate = rep(s, each = 8), #add PCR plate information
         RNAplate = rep(1:length(s), each = 8), #add RNA plate information
         row = rep(1:8, nrow(data) / 8), #add to which PCR-plate the rows belong
         group = rep(1:(nrow(data) / 16), each = 16)) %>%  #add a row number
  gather(key = "column", value = "wells", X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12) %>% #gather all columns 
  mutate(id = 1:n()) #id all wells

all_data <- d %>%
  full_join(dplyr::select(w, id, wells), by = c("id" = "id")) #combine Cq and well information


#Define filters for dataset
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

avg_all <- mean(all_data %>% #calculate the mean of all control values (filtered for full plates with low prevalence and good controls)
                  filter(!(RNAplate %in% many_plate), !(RNAplate %in% nonfull_plate),  !(RNAplate %in% bad_control), row == 8, column == "X12") %>%
                  .$cq, na.rm = TRUE)

correction <- all_data %>% 
  filter(row == 8 & column == "X12") %>%
  group_by(PCRplate) %>% #calculate the mean of the positive controls on each PCR-plate
  summarize(avg = mean(cq, na.rm = TRUE), diff = avg - avg_all) %>% #find how much the values in this plate should be adjusted
  ungroup()


min = mean(correction$avg, na.rm = TRUE) - 2*sd(correction$avg, na.rm = TRUE) #calculate the lower boundary
max = mean(correction$avg, na.rm = TRUE) + 2*sd(correction$avg, na.rm = TRUE) #calculate the upper boundary

good_cont <- correction %>% #Find plates whose positive control value is between min and max
  filter(avg > min, avg < max) %>%
  .$PCRplate


#Filter out all plates with have too many positives, are not full or have bad control values. Also remove any control wels
d_fil <- all_data %>%
  filter(!(RNAplate %in% many_plate), !(id %in% con_ids), !(RNAplate %in% nonfull_plate),  !(RNAplate %in% bad_control), PCRplate %in% good_cont) %>%
  left_join(correction, by = c("PCRplate" = "PCRplate")) %>% #join the correction
  mutate(corr_cq = cq - diff) %>% #calculate the correction
  arrange(PCRplate, RNAplate, row) #just to make it nicer


#Get the corrected cq values. Since 37 is the cut-off, we will only use values under 37
cq_pos <- d_fil %>%
  filter(corr_cq <= 37) %>%
  .$corr_cq

#Data wrangling for RDML publication
data_pub <- d_fil %>%
  filter(corr_cq <= 37)  %>%
  mutate(target = "E gene", sample_type = "std", plate_id = PCRplate,
         target_type = "toi", sample = id, dye = "FAM", reaction = ifelse(plate_id %in% plates$X1[1:155], "singleplex", "duplex")) %>%
  dplyr::select(plate_id, sample, sample_type, target, target_type, dye, reaction, wells, corr_cq)
write_csv(data_pub, "../../Data/filtered_cq.csv")

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

#Function to call for simulation
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

######### Simulation ############

#Define strategies
pools <- list(c(1,4), c(1,8), c(1, 12), c(1, 16), c(1, 24), c(8, 12), c(12, 16), c(16, 24))

#Run simulations
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

#Combine everything
s <- rbind(s1 %>% mutate(rep = 1), 
           s2 %>% mutate(rep = 2),
           s3 %>% mutate(rep = 3), 
           s4 %>% mutate(rep = 4),
           s5 %>% mutate(rep = 5)) %>%
  mutate(sensitivity = 1 - falseneg)

#Define the names of the strategies
strategies <- c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")

#Collect all Cq values from the positive samples
p <- rbind(list.rbind(s1$positive_values),
                  list.rbind(s2$positive_values),
                  list.rbind(s3$positive_values),
                  list.rbind(s4$positive_values),
                  list.rbind(s5$positive_values)) %>%
  mutate(rep = rep(c("1", "2", "3", "4", "5"), each = nrow(.) / 5),
         strategy = factor(rep(rep(strategies, each = 91871), 5), levels = strategies),
         prevalence = test / 1000,
         tested = as.factor(ifelse(tested, "TP", "FN"))) 


############ Figures ###############
#Figure 1
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

#Figure 2
s %>%
  mutate(efficiency = 100000 / ntest, 
         group = factor(paste0(row, "x", col), levels = strategies),
         Prevalence = pos / 100000 * 100) %>%
  group_by(group, Prevalence) %>%
  summarise(sens = median(sensitivity), eff = median(efficiency)) %>%
  ungroup() %>%
  ggplot(aes(Prevalence, sens, col = eff)) +
  geom_point() +
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

#Figure 3
p %>%
  filter(prevalence %in% c(0.01, 0.10, 1, 10)) %>%
  group_by(strategy, prevalence, posneg.x) %>%
  summarise(TP = sum(tested == "TP"), total = n()) %>%
  ungroup() %>%
  group_by(strategy, prevalence) %>%
  mutate(cs = cumsum(TP), cs_total = cumsum(total), sens = cs / cs_total) %>%
  ggplot(aes(posneg.x, sens, col = strategy)) +
  geom_line() +
  facet_wrap(~prevalence, 
             strip.position = "bottom",
             nrow = 1) +
  theme_minimal() +
  labs(caption = "prevalence (%)") +
  ylab("cumulative sensitivity\n") +
  xlab("Cq\n") +
  theme(axis.title = element_text(face = "bold", size  = 15),
        strip.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"),
        plot.caption = element_text(face = "bold", size = 15, hjust = 0.5)) +
  scale_colour_manual(values=cbPalette) +
  scale_x_continuous(position = "top", limits = c(30, 37)) 

#Supplemental Figure 1
d_fil %>%
  mutate(group = ifelse(corr_cq > 37, "grey", "white")) %>%
  ggplot(aes(corr_cq, fill = group)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), col = "#2D5580", binwidth = 1, size = 1, breaks = 9:42) +
  xlab("\nCq value") +
  ylab("fraction\n") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(10, 40, by = 5)) +
  scale_fill_manual(values = c("grey", "white"))

#Supplemental Figure 2
p %>%
  ggplot(aes(tested, posneg.x, col = strategy)) +
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

#Supplemental Figure 3
s %>%
  mutate(efficiency = 100000 / ntest, group = factor(paste0(row, "x", col), levels = c("1x4", "1x8", "1x12", "1x16", "1x24", "8x12", "12x16", "16x24")),
         Prevalence = pos / 100000 * 100,
         dimension = ifelse(row == 1, "1D", "2D")) %>%
  group_by(group, Prevalence) %>%
  summarise(sens = median(sensitivity)) %>%
  ungroup() %>%
  ggplot(aes(group, sens)) +
  geom_boxplot(col = "black") +
  geom_beeswarm(aes(col = Prevalence)) +
  ylab("sensitivity\n") +
  xlab("") +
  ylim(0.5, 1.02) +
  labs(col = "prevalence\n") +
  theme_minimal() +
  scale_color_gradientn(colours = gradient,
                        trans = "log", 
                        breaks = c(0.01, 0.1, 1, 10)) +
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
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.margin=unit(c(0,4.5,0,1),"cm")) 

#Supplemental Figure 4
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
             nrow = 2,
             strip.position = "bottom",
             scales = "free_x") +
  theme_minimal() +
  ylab("cumulative sensitivity\n") +
  xlab("Cq\n") +
  xlim(30, 37) +
  theme(axis.title = element_text(face = "bold", size  = 15),
        strip.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  scale_colour_manual(values=cbPalette[]) +
  scale_x_continuous(position = "top", limits = c(30, 37)) 

