# PREAMBLE ----

## LIBRARIES ----

### Data Management ----
require(argparse)
require(tidyverse)
require(dplyr)
require(data.table)
require(zoo)

### Visualization and Plotting ----
require(ggpubr)
require(ggtext)
require(ggrepel)
require(scales)
require(cowplot)
require(corrplot)
require(colorspace)
require(viridis)
require(PerformanceAnalytics)
require(wesanderson)
require(jtools)

### Dates and Times ----
require(anytime)
require(lubridate)

### Free Text and Strings ----
require(stringr)
require(stringi)
require(stringdist)
require(english)
require(tm)
require(topicmodels)
require(lattice)
require(tidytext)
require(proxy)
require(corpus)
require(gender)
require(rxnorm)
require(RxNormR)
require(qdap)
data(stop_words)

### Locations ----
require(zipcodeR)

### Correlation and Regression ----
require(pspearman)
require(jtools)
require(caret)

### Missing Data ----
require(naniar)
require(finalfit)

### Dimensional Reduction ----
require(psych)
require(nFactors)
require(FactoMineR)
require(factoextra)
require(corrr)
require(car)
require(BBmisc)
require(mice)
require(phateR)


## ARGUMENTS ----
parser <- ArgumentParser()
parser$add_argument("--dir",
                    help = "Working directory",
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data")

parser$add_argument("--date",
                    help = "Data freeze date",
                    default = "20221120")

parser$add_argument("--geno",
                    help = "Genetic dataset IDs",
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/gwa/gen/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465.fam")

args <- parser$parse_args()
print(args)

### Set working directory ----
workDir = args$dir
setwd(workDir)

### Set data freeze ----
freezeDate = args$date
freezeDir = paste(workDir, "data", "database", "data_freeze", freezeDate, sep = "/")

### Load from saved image, if available ----
load(paste(workDir,
           "dat",
           paste(
             paste("DarwinsArk",
                   freezeDate,
                   "phenotyping", sep = "_"),
             ".RData",
             sep = ""
           ), sep = "/"))

## LOAD DATA ----
dogs = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_dogs.csv",
  sep = ""
))

questions = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_questions.csv",
  sep = ""
))

answers = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_answers.csv",
  sep = ""
))

scores = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_factor-scores.csv",
  sep = ""
))

gwas = read_tsv(args$geno,
                col_names = c("fam",
                              "dog",
                              "pat",
                              "mat",
                              "sex",
                              "phe"))

anc = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/dat/adm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_global-adm_supervised_K-115.csv")

# PHENOTYPE CODES ----

# Scores from factor analysis
# fa.{code}

# Responses to survey items
# bq.{code}

# Morphology phenotypes (discrete)
# mp.{code}

# Morphology phenotypes (quantitative)
# mq.{code}

# Scores from allergy PCA
# af.{code}

# TYPICAL PHENOTYPES ----

# Phenotype domains:
# fa - factor analysis
# bq - survey item
# mp - morphology, discrete
# mq - morphology, quantitative
# af - allergy factor

## FACTOR SCORES: ----

# `fa.[ID].filled-by-mean`: quantitative phenotype, factor scores as a linear combination of scaled survey responses from dogs with >=75% completeness of component items
# `fa.[ID].tail`: binary phenotype, cases as lower tail (Q25) of factor scores and controls as all else

scores = scores %>%
  mutate(factor = paste(
    "fa.",
    str_pad(
      string = as.character(factor),
      width = 2,
      side = "left",
      pad = "0"
    ),
    sep = ""
  ))

### Extract phenotypes and covariates ----

fa.ids = unique(scores$factor)

fa.all=list()
i=1
for (id in fa.ids){
  
  # get scores
  scores.id = scores %>%
    filter(factor == id) %>%
    select(dog,
           age=age_mean,
           phe=score_fill_mean_norm) %>%
    merge((dogs %>%
             select(dog,
                    hgt=height_filled)), 
          by = "dog", 
          all.x = T) %>%
    mutate(fid=dog,
           iid=dog) %>%
    select(fid,
           iid,
           phe,
           age,
           hgt,
           dog)
  
  # extract age covariates:
  write_tsv(x = (scores.id %>%
                   filter(dog %in% gwas$dog) %>%
                   select(fid,iid,age)),
            paste(workDir,
                  "/gwa/cov/qcov/",
                  "age.",
                  id,
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract age + hgt covariates:
  write_tsv(x = (scores.id %>%
                   filter(dog %in% gwas$dog) %>%
                   select(fid,iid,age,hgt)),
            paste(workDir,
                  "/gwa/cov/qcov/",
                  "age.hgt.",
                  id,
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract fa.[ID]-filled-by-mean phenotypes:
  fa.phe = scores.id %>%
    filter(dog %in% gwas$dog) %>%
    select(fid,iid,phe) %>%
    mutate(phe = round(phe,
                       digits = 2))
  
  fa.all[[i]] = fa.phe %>%
    mutate(code = paste(id,
                        "-filled-by-mean",
                        sep = ""))
  
  write_tsv(x = fa.phe,
            paste(workDir,
                  "/gwa/phe/",
                  id,
                  "-filled-by-mean",
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract `fa.[ID].tail` phenotypes:
  fa.phe =scores.id %>%
    mutate(Q25 = quantile(phe,0.25)) %>%
    mutate(phe = if_else(phe<Q25,
                         1,
                         0)) %>%
    filter(dog %in% gwas$dog) %>%
    select(fid,iid,phe)
  
  fa.all[[i]] = fa.all[[i]] %>%
    bind_rows(fa.phe %>%
    mutate(code = paste(id,
                        ".tail",
                        sep = "")))
  
  write_tsv(x = fa.phe,
            paste(workDir,
                  "/gwa/phe/",
                  id,
                  ".tail",
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract `fa.[ID].tail.fill` phenotypes:
  fa.phe = scores.id %>%
    mutate(Q25 = quantile(phe,0.25)) %>%
    mutate(phe = if_else(phe<Q25,
                         1,
                         0)) %>%
    filter(dog %in% gwas$dog) %>%
    select(fid,iid,phe)
  
  fa.phe = bind_rows(fa.phe,
                     data.table(fid= (gwas %>%
                                        filter(!dog %in% fa.phe$iid) %>%
                                        pull(dog)),
                                iid = (gwas %>%
                                         filter(!dog %in% fa.phe$iid) %>%
                                         pull(dog)),
                                phe=0)) %>%
    arrange(iid)
  
  fa.all[[i]] = fa.all[[i]] %>%
    bind_rows(fa.phe %>%
                mutate(code = paste(id,
                                    ".tail.fill",
                                    sep = "")))
  
  write_tsv(x = fa.phe,
            paste(workDir,
                  "/gwa/phe/",
                  id,
                  ".tail.fill",
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  i=i+1
}

### Generate `fa.all` table ----
fa.all = rbindlist(fa.all) %>%
  pivot_wider(id_cols = c("fid","iid"),
              names_from = "code",
              values_from = "phe")

### Plot distributions ----
for (fa in unique(scores$factor)) {
  p = bind_rows((scores %>% 
                   filter(factor == fa) %>%
                   mutate(Q25 = quantile(score_fill_mean_norm, 0.25)) %>%
                   mutate(phe = round(score_fill_mean_norm, 
                                      digits = 2)) %>%
                   filter(!is.na(phe)) %>%
                   mutate(n = n()) %>%
                   mutate(set = paste("all dogs ",
                                      "(n= ",n,")",
                                      sep = "")) %>%
                   select(dog,phe,set,Q25)),
                (scores %>% 
                   filter(factor == fa) %>%
                   mutate(Q25 = quantile(score_fill_mean_norm, 0.25)) %>%
                   mutate(phe = round(score_fill_mean_norm, 
                                      digits = 2)) %>%
                   filter(!is.na(phe)) %>%
                   filter(dog %in% gwas$dog) %>%
                   mutate(n = n()) %>%
                   mutate(set = paste("genetic data ",
                                      "(n= ",n,")",
                                      sep = "")) %>%
                   select(dog,phe,set,Q25))) %>%
    ggplot(aes(x = phe)) +
    geom_histogram(aes(fill = set, y=..density..),
                   alpha = 0.25, binwidth = 0.25,
                   position="identity") +
    geom_density(aes(color = set, linetype = set),
                 linewidth = 0.75) +
    geom_vline(aes(xintercept = Q25),
               linewidth = 2, alpha = 0.25) +
    labs(x = "score",
         caption = paste("phenocode:",fa)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_discrete(type = rev(greenfocus[1:2])) +
    scale_fill_discrete(type = lighten(rev(greenfocus[1:2]),amount = .4)) +
    scale_linetype_discrete() +
    theme_pubr()
  
  ggsave(filename = paste(workDir,
                          "/gwa/phe/",
                          "dist-plot.",
                          fa,
                          ".png",
                          sep = ""),
         plot = p,
         device = "png",
         dpi = 300,
         units = "in",
         width = 3*2,
         height = 3*2)
  
  ggsave(filename = paste(workDir,
                          "/gwa/phe/",
                          "dist-plot.",
                          fa,
                          ".pdf",
                          sep = ""),
         plot = p,
         device = "pdf",
         dpi = 300,
         units = "in",
         width = 3*2,
         height = 3*2)
}

### Breed-aggregate sampled scores ----
fa.summ = scores %>%
  merge((dogs %>%
           select(dog,breed)), 
        by = "dog", 
        all.x = T) %>%
  mutate(score = score_fill_mean) %>%
  group_by(breed,factor) %>%
  summarize(score.mean = mean(score, na.rm = T),
            score.sd = sd(score, na.rm = T))

fa.summ.all = scores %>%
  group_by(factor) %>%
  mutate(score = score_fill_mean) %>%
  summarize(score.mean = mean(score, na.rm = T),
            score.sd = sd(score, na.rm = T))

fa.list = unique(fa.summ$factor)
fa.samp.list = list()
i=1
for (fa in fa.list){
  fa.phe = anc %>%
    group_by(dog) %>%
    slice_max(pct, n=1) %>%
    mutate(purebred = pct > 0.85) %>%
    mutate(breed = if_else(purebred,
                           pop,
                           NA_character_)) %>%
    mutate(score = if_else(is.na(breed),
                           rnorm(mean=(fa.summ.all %>% filter(factor == fa) %>% pull(score.mean)),
                                 sd=(fa.summ.all %>% filter(factor == fa) %>% pull(score.sd)),
                                 n=1),
                           rnorm(mean=(fa.summ %>% filter(breed==breed & factor == fa) %>% pull(score.mean)),
                                 sd=(fa.summ %>% filter(breed==breed & factor == fa) %>% pull(score.sd)),
                                 n=1))) %>% 
    ungroup() %>%
    mutate(score = round(score, digits = 2)) 

  fa.samp.list[[i]] = fa.phe %>%
    mutate(factor=paste("fa.breed-agg.",
                        fa,
                        sep = "")) %>%
    select(dog,factor,score)

  write_tsv(x = (fa.phe %>%
                   mutate(fid=dog,
                          iid=dog,
                          phe=score) %>%
                   select(fid,iid,phe)),
            paste(workDir,
                  "/gwa/phe/",
                  "fa.breed-agg.",
                  fa,
                  ".tsv",
                  sep = ""),
            col_names = F)

  i=i+1
}

fa.agg.all = rbindlist(fa.samp.list) %>%
  pivot_wider(id_cols = "dog",
              names_from = "factor",
              values_from = "score") %>%
  mutate_all(~ifelse(is.nan(.), NA, .))

## SURVEY ITEMS: ----
# `bq.[ID]`: quantitative phenotype, responses to ordinal survey items as-is
# `bq.[ID].mean-binary`: binary phenotypes,

bq.data = answers %>%
  filter(question %in%
           c(1:110,
           144:155,
           165:175,
           177:182,
           184,185,
           195:240)) %>%
  select(question,dog,answer,option,age) %>%
  mutate(answer = if_else(option %in% c("I don't know",
                                        "I'm not sure",
                                        "Maybe"),
                          "",
                          answer)) %>%
  mutate(trait = paste("bq.",
                       str_pad(question,
                               3,
                               "left",
                               "0"), 
                       sep = "")) %>%
  mutate(fid = dog,
         iid = dog,
         phe = as.numeric(answer)) %>%
  merge((dogs %>%
           select(dog,
                  hgt=height_filled)),
        by = "dog", all.x = T) %>%
  select(trait,fid,iid,phe,age,hgt,dog)

### Extract phenotypes ----
bq.ids = unique(bq.data$trait)

bq.all=list()
i=1

for (id in bq.ids) {
  
  # get data
  bq.data.id = bq.data %>%
    filter(trait == id) %>%
    select(fid,
           iid,
           phe,
           age,
           hgt,
           dog)
  
  # extract age covariates:
  
  # extract age covariates:
  write_tsv(x = (bq.data.id %>%
                   filter(dog %in% gwas$dog) %>%
                   select(fid,iid,age)),
            paste(workDir,
                  "/gwa/cov/qcov/",
                  "age.",
                  id,
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract age + hgt covariates:
  write_tsv(x = (bq.data.id %>%
                   filter(dog %in% gwas$dog) %>%
                   select(fid,iid,age,hgt)),
            paste(workDir,
                  "/gwa/cov/qcov/",
                  "age.hgt.",
                  id,
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract `bq.[ID]` phenotypes:
  bq.phe = bq.data.id %>%
    filter(dog %in% gwas$dog) %>%
    select(fid,iid,phe)
  
  bq.all[[i]] = bq.phe %>%
    mutate(code = id)
  
  write_tsv(x = bq.phe,
            paste(workDir,
                  "/gwa/phe/",
                  id,
                  ".tsv",
                  sep = ""),
            col_names = F)

  # extract `bq.[ID].mean-binary` phenotypes:
  bq.phe = bq.data.id %>%
    mutate(mean = mean(phe, na.rm = T)) %>%
    mutate(phe = if_else(phe > mean,
                         1,
                         0)) %>%
    filter(dog %in% gwas$dog) %>%
    select(fid,iid,phe)
  
  bq.all[[i]] = bind_rows(bq.all[[i]],
                          bq.phe %>%
    mutate(code = paste(id,
                        ".mean-binary",
                        sep = "")))
  
  write_tsv(x = bq.phe,
            paste(workDir,
                  "/gwa/phe/",
                  id,
                  ".mean-binary",
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  # extract `bq.[ID].mean.binary.fill` phenotypes:
  bq.phe = bq.data.id %>%
    mutate(mean = mean(phe, na.rm= T)) %>%
    mutate(phe = if_else(phe > mean,
                         1,
                         0)) %>%
    filter(dog %in% gwas$dog) %>%
    select(fid,iid,phe)
  
  bq.phe = bind_rows(bq.phe,
                     data.table(fid= (gwas %>%
                                        filter(!dog %in% bq.phe$iid) %>%
                                        pull(dog)),
                                iid = (gwas %>%
                                         filter(!dog %in% bq.phe$iid) %>%
                                         pull(dog)),
                                phe=0)) %>%
    arrange(iid)
  
  write_tsv(x = bq.phe,
            paste(workDir,
                  "/gwa/phe/",
                  id,
                  ".tail.binary.fill",
                  ".tsv",
                  sep = ""),
            col_names = F)
  
  i=i+1
}

### Generate `bq.all` table ----
bq.all = rbindlist(bq.all) %>%
  pivot_wider(id_cols = c("fid","iid"),
              names_from = "code",
              values_from = "phe")

### Plot distributions ----
for (q in c(1:110,
             144:155,
             165:175,
             177:182,
             184,185,
             195:240)) {
  bq = paste("bq.",str_pad(q,3,"left","0"), sep = "")
  
  thisQuestion = answers %>%
    filter(question == q) %>%
    select(question,dog,answer,option,age) %>%
    mutate(answer = if_else(option %in% c("I don't know",
                                          "I'm not sure",
                                          "Maybe"),
                            "",
                            answer)) %>%
    mutate(fid = dog,
           iid = dog,
           phe = as.numeric(answer)) %>%
    merge((dogs %>%
             select(dog,
                    hgt=height_filled)),
          by = "dog", all.x = T)
  
  p = bind_rows((thisQuestion %>%
                   filter(!is.na(phe)) %>%
                   mutate(n = n()) %>%
                   mutate(set = paste("all dogs ",
                                      "(n= ",n,")",
                                      sep = "")) %>%
                   select(dog,phe,set,option)),
                (thisQuestion %>%
                   filter(!is.na(phe)) %>%
                   filter(dog %in% gwas$dog) %>%
                   mutate(n = n()) %>%
                   mutate(set = paste("genetic data ",
                                      "(n= ",n,")",
                                      sep = "")) %>%
                   select(dog,phe,set,option))) %>%
    mutate(option = if_else(option == "Neither Agree Nor Disagree",
                            "Neither",
                            option)) %>%
    mutate(lab = paste(as.character(as.integer(phe)),
                       ". ", option,
                       sep = "")) %>%
    ungroup() %>%
    mutate(mean_phe = mean(as.integer(phe), na.rm = T)) %>%
    group_by(phe,lab,set,mean_phe) %>%
    summarize(phe_n = n()) %>%
    mutate(mean_opt = if_else(phe == round(mean_phe),
                              lab,
                              NA_character_)) %>%
    ggplot(aes(x = lab,
               y = phe_n)) +
    coord_flip(expand = F, 
               clip = F) +
    facet_wrap(.~set,
               nrow = 2) +
    geom_vline(aes(xintercept = mean_opt),
               linewidth = 2, alpha = 0.25) +
    geom_col(aes(color = set, fill = set),
                 linewidth = 0.75, width = 0.5) +
    geom_text(aes(label = phe_n),
              hjust = -0.25,
              size = 2,
              position = position_dodge(width = 1.25)) +
    labs(x = "option",
         y = "# dogs",
         caption = paste("phenocode:",bq),
         title = (questions %>% filter(id == q) %>% pull(string))) +
    scale_y_continuous(limits = c(0,(thisQuestion %>%
                                       group_by(option) %>%
                                       summarize(n = n()) %>%
                                       slice_max(order_by = n, n = 1) %>% pull(n))+500)) +
    scale_color_discrete(type = rev(greenfocus[1:2])) +
    scale_fill_discrete(type = lighten(rev(greenfocus[1:2]),amount = .4)) +
    scale_linetype_discrete() +
    theme_pubr() +
    theme(legend.position="bottom", 
          axis.title.y = element_blank(),
          plot.margin = margin(0,1,0,0, "cm"), plot.title = element_text(size=10))
  
  ggsave(filename = paste(workDir,
                          "/gwa/phe/",
                          "dist-plot.",
                          bq,
                          ".png",
                          sep = ""),
         plot = p,
         device = "png",
         dpi = 300,
         units = "in",
         width = 3*2,
         height = 3*2)
  
  ggsave(filename = paste(workDir,
                          "/gwa/phe/",
                          "dist-plot.",
                          bq,
                          ".pdf",
                          sep = ""),
         plot = p,
         device = "pdf",
         dpi = 300,
         units = "in",
         width = 3*2,
         height = 3*2)
}

## MORPHOLOGY: ----
# get weight separate:
mq.weight = answers %>%
  filter(question == 242) %>%
  select(dog,answer,age) %>%
  mutate(answer = gsub("-lb"," lb", gsub("-kg", " kg", answer))) %>%
  separate(answer, into = c('value','unit'), sep = " ", convert = T) %>%
  mutate(value = abs(value)) %>%
  mutate(weight = if_else(unit == "lb",
                          value,
                          value*2.2)) %>%
  select(dog, weight)

# if multichoice, then split:
mq.multi = answers %>%
  as_tibble() %>%
  filter(question %in% c(121:128,243:250)) %>%
  select(question,dog,option) %>%
  filter(question %in% (questions %>% 
                          filter(style %in% c("MultiChoices",
                                              "MultiChoicesWithoutNone",
                                              "MultiChoicesTooltip",
                                              "Colors")) %>% 
                          pull(id))) %>%
  filter(!option %in% c("I don't know",
                        "I'm not sure",
                        "Maybe")) %>%
  mutate(option = gsub(" ",
                       "_",
                       tolower(str_trim(gsub("\\s*\\([^\\)]+\\)",
                                             "",
                                             option))))) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = c("question","option"),
              names_sep = "_",
              values_from = "value",
              values_fn = min) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate_if(is.integer, ~replace_na(., 0))

# else, keep
mq.quant = answers %>%
  as_tibble() %>%
  filter(question %in% c(121:128,243:250)) %>%
  select(question,dog,option,answer) %>%
  filter(!question %in% (questions %>% filter(style %in% c("MultiChoices",
                                                           "MultiChoicesWithoutNone",
                                                           "MultiChoicesTooltip",
                                                           "Colors")) %>% pull(id))) %>%
  filter(!option %in% c("I don't know",
                        "I'm not sure",
                        "Maybe")) %>%
  select(dog,question,answer) %>%
  unique() %>%
  mutate(answer = as.numeric(answer)) %>%
  pivot_wider(names_from = "question",
              names_prefix = "Q",
              values_from = "answer",
              values_fn = min) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate_if(is.integer, ~replace_na(., 0))

mq = mq.weight %>% 
  merge(mq.quant, by = "dog", all = T) %>% 
  merge(mq.multi, by = "dog", all = T) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate_if(is.integer, ~replace_na(., 0))

# morphology questions (mq):
mq.125 = mq.quant$Q125
mq.125[mq.125 %in% c(5,6)] = NA
mq.quant$Q125 = mq.125

mq.all = mq.weight %>% 
  merge(mq.quant, by = "dog", all = T) %>%
  rename(Q242=weight) %>%
  filter(dog %in% gwas$dog)
colnames(mq.all) = gsub("Q","mq.",colnames(mq.all))

for (P in colnames(mq.all)[-1]){
  write_tsv(x = (mq.all %>%
                   filter(dog %in% gwas$dog) %>%
                   select(dog,!!P) %>%
                   mutate(fid=dog,
                          iid=dog,
                          phe=get(P)) %>%
                   filter(!is.na(phe)) %>%
                   select(fid,iid,phe) %>%
                   arrange(iid) %>%
                   as.data.table()),
            paste(workDir,
                  "/gwa/phe/",
                  P,
                  ".tsv",
                  sep = ""),
            col_names = F)
}

# eye color
mp_eye_color = answers %>%
  filter(question == 249) %>%
  as.data.table()

mp_eye_color[answer == 4, phe := 3]
mp_eye_color[answer == 3, phe := 3]
mp_eye_color[answer == 0, phe := 2]
mp_eye_color[answer == 2, phe := 1]
mp_eye_color[answer == 1, phe := 0]

mp_eye_color = mp_eye_color %>% 
  select(dog,`mp.249.eye-color`=phe) %>%
  group_by(dog) %>%
  slice_head(n=1) %>%
  ungroup()

# odd eyes
mp_eye_heter = bind_rows((answers %>%
                            filter(question == 249) %>%
                            select(dog,option) %>%
                            group_by(dog) %>%
                            mutate(n = length(unique(option))) %>%
                            group_by(dog) %>%
                            summarize(phe = n > 1) %>%
                            ungroup() %>%
                            unique()),
                         (answers %>%
                            filter(question == 126) %>%
                            select(dog,option) %>%
                            filter(option != "I don't know") %>%
                            mutate(phe = if_else(option == "Yes",
                                   1,
                                   0)) %>%
                            select(dog,phe) %>%
                            unique())) %>%
  select(dog,`mp.126.249.eye-heter`=phe) %>%
  ungroup() %>%
  group_by(dog) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  unique()

# pad color
mp_pad_color = answers %>%
  filter(question == 250) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.250.pad-liver` = `1`,
         `mp.250.pad-black` = `0`,
         `mp.250.pad-grey` = `2`,
         `mp.250.pad-pink` = `3`)

# features
mp_features = answers %>%
  filter(question == 248) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.248.jowled` = `0`,
         `mp.248.wrinkly` = `1`,
         `mp.248.freckled` = `2`,
         `mp.248.full_white` = `3`,
         `mp.248.vitiligo` = `4`,
         `mp.248.ridged` = `5`,
         `mp.248.hairless` = `6`,
         `mp.248.underbite` = `7`,
         `mp.248.tongue_color` = `8`) %>%
  select(-`9`)

# head shape
mp_head_shape = answers %>%
  filter(question == 245) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  mutate(`mp.245.brachycephaly` = if_else(`2`==1 | `4`==1 | `5`==1,
                                          1,
                                          0),
         `mp.245.dolichocephaly` = if_else(`0`==1 | `1`==1,
                                           1,
                                           0)) %>%
  rename(`mp.245.head-shape-A` = `0`,
         `mp.245.head-shape-B` = `1`,
         `mp.245.head-shape-C` = `2`,
         `mp.245.head-shape-D` = `3`,
         `mp.245.head-shape-E` = `4`,
         `mp.245.head-shape-F` = `5`)

mp_coat_colors_243 = answers %>%
  filter(question == 243) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.243.coat-color-black` = `0`,
         `mp.243.coat-color-brown` = `1`,
         `mp.243.coat-color-white` = `2`,
         `mp.243.coat-color-red` = `3`,
         `mp.243.coat-color-yellow` = `4`,
         `mp.243.coat-color-dilute-black` = `5`,
         `mp.243.coat-color-dilute-brown` = `6`,
         `mp.243.coat-color-skin` = `7`,
         `mp.243.coat-color-tan` = `8`,
         `mp.243.coat-color-cream` = `9`)

mp_coat_colors_122 = answers %>%
  filter(question == 122) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.122.coat-color-white` = `0`,
         `mp.122.coat-color-red` = `1`,
         `mp.122.coat-color-yellow` = `2`,
         `mp.122.coat-color-gray` = `3`,
         `mp.122.coat-color-chocolate-brown` = `4`,
         `mp.122.coat-color-pure-black` = `5`,
         `mp.122.coat-pattern-merle` = `6`,
         `mp.122.coat-pattern-brindle` = `7`) %>%
  select(-`8`)

mp_coat_patterns = answers %>%
  filter(question == 244) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.244.coat-pattern-striping` = `0`,
         `mp.244.coat-pattern-patching` = `1`,
         `mp.244.coat-pattern-splotching` = `2`,
         `mp.244.coat-pattern-spotting` = `3`,
         `mp.244.coat-pattern-roaning` = `4`,
         `mp.244.coat-pattern-ticking` = `5`,
         `mp.244.coat-pattern-none` = `6`)

mp.coat_patterns_spots = answers %>%
  filter(question == 244) %>%
  select(dog,answer) %>%
  filter(answer %in% c(3,4,5)) %>%
  select(dog,answer) %>%
  as.data.table()
mp.coat_patterns_spots[answer==3, `mq.244.coat-pattern-spots` := 3]
mp.coat_patterns_spots[answer==4, `mq.244.coat-pattern-spots` := 2]
mp.coat_patterns_spots[answer==5, `mq.244.coat-pattern-spots` := 1]

mp.coat_patterns_spots = mp.coat_patterns_spots %>%
  group_by(dog) %>%
  summarize(`mq.244.coat-pattern-spots` = max(`mq.244.coat-pattern-spots`)) %>%
  bind_rows((answers %>%
               filter(question == 244) %>%
               filter(!dog %in% mp.coat_patterns_spots$dog) %>%
               select(dog) %>%
               unique() %>%
               mutate(`mq.244.coat-pattern-spots` = 0))) %>%
  arrange(dog)

mp.earshape = answers %>%
  filter(question == 125) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.125.ears-pendant` = `0`,
         `mp.125.ears-dropped` = `1`,
         `mp.125.ears-rose-button` = `2`,
         `mp.125.ears-pricked` = `3`,
         `mp.125.ears-cropped` = `4`,
         `mp.125.ears-other` = `5`)

mp.all = mp_eye_color %>%
  merge(mp_eye_heter, by = "dog", all = T) %>%
  merge(mp_features,by = "dog", all = T) %>%
  merge(mp_head_shape,by = "dog", all = T) %>%
  merge(mp.earshape,by = "dog", all = T) %>%
  merge(mp_pad_color,by = "dog", all = T) %>%
  merge(mp_coat_colors_243,by = "dog", all = T) %>%
  merge(mp_coat_colors_122,by = "dog", all = T) %>%
  merge(mp_coat_patterns,by = "dog", all = T) %>%
  merge(mp.coat_patterns_spots,by = "dog", all = T) %>%
  filter(dog %in% gwas$dog) %>%
  unique()

for (P in colnames(mp.all)[-1]){
  write_tsv(x = (mp.all %>%
                   select(dog,!!P) %>%
                   mutate(fid=dog,
                          iid=dog,
                          phe=get(P)) %>%
                   filter(!is.na(phe)) %>%
                   select(fid,iid,phe) %>%
                   arrange(iid) %>%
                   as.data.table() %>%
                   unique()),
            paste(workDir,
                  "/gwa/phe/",
                  P,
                  ".tsv",
                  sep = ""),
            col_names = F)
}

# ATYPICAL PHENOTYPES ----

## LINEAGE PSEUDOTIME VALUES for 12 lineages ----
lineage.ext = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/gwa/phe/lineage.ext.tsv")
colnames(lineage.ext) = gsub(" ","",colnames(lineage.ext))

for (P in colnames(lineage.ext)[-1]){
  write_tsv(x = (lineage.ext %>%
                   select(dog,!!P) %>%
                   mutate(fid=dog,
                          iid=dog,
                          phe=get(P)) %>%
                   filter(!is.na(phe)) %>%
                   select(fid,iid,phe) %>%
                   arrange(iid) %>%
                   as.data.table()),
            paste(workDir,
                  "/gwa/phe/",
                  P,
                  ".tsv",
                  sep = ""),
            col_names = F)
}

lineage.int = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/gwa/phe/lineage.int.tsv")
colnames(lineage.int) = gsub(" ","",colnames(lineage.int))


for (P in colnames(lineage.int)[-1]){
  write_tsv(x = (lineage.int %>%
                   select(dog,!!P) %>%
                   mutate(fid=dog,
                          iid=dog,
                          phe=get(P)) %>%
                   filter(!is.na(phe)) %>%
                   select(fid,iid,phe) %>%
                   arrange(iid) %>%
                   as.data.table()),
            paste(workDir,
                  "/gwa/phe/",
                  P,
                  ".tsv",
                  sep = ""),
            col_names = F)
}


## BREED-AGGREGATE C-BARQ SCORES from 12,000 dogs of 203 breeds ----
cbarq = read_tsv("~/Dropbox (UMass Medical School)/Thesis/data/CBARQ.dat.tsv") %>%
  mutate(breed = gsub(" ","_",tolower(BreedID)))

cb.list = cbarq %>% 
  select(-c("BreedID","breed")) %>%
  colnames()

### Calculated per-breed mean and sd ----
cb.summ = cbarq %>%
  select(-BreedID) %>%
  pivot_longer(cols = !breed,
               names_to = "factor",
               values_to = "score") %>%
  group_by(breed,factor) %>%
  summarize(score.mean = mean(score, na.rm = T),
            score.sd = sd(score, na.rm = T))

### Calculate breed-wide mean and sd ----

# or don't re-run, just load:
# cb.all=list()
# i=1
# for (cb in unique(cb.summ$factor)){
#   cb.all[[i]] = read_tsv(paste(workDir,
#                                "/gwa/phe/",
#                                "cbarq.breed-sampled.",
#                                cb,
#                                ".tsv",
#                                sep = ""),
#                          col_names = c("fid","iid","phe"))
#   cb.all[[i]]$factor = paste("cbarq.breed-sampled.",
#                              cb,
#                              sep = "")
#   i=i+1
# }

# cb.all = rbindlist(cb.all) %>%
#   pivot_wider(id_cols = c("fid","iid"),
#               names_from = "factor",
#               values_from = "phe")

cb.all = cbarq %>%
  select(-c(BreedID,breed)) %>%
  pivot_longer(cols = everything(),
               names_to = "factor",
               values_to = "score") %>%
  group_by(factor) %>%
  summarize(score.mean = mean(score, na.rm = T),
            score.sd = sd(score, na.rm = T))

set.seed(1337435555)

### Plot normal C-BARQ distributions ----
cb.norms = cb.summ %>%
  filter(breed %in% (anc %>% filter(pct > 0.85) %>% group_by(pop) %>% summarize(n = n()) %>% filter(n>10) %>% pull(pop))) %>%
  group_by(breed,factor) %>%
  summarize(val = rnorm(10000,
                        mean = score.mean,
                        sd = score.sd)) %>%
  ungroup()

p = cb.norms %>%
  ggplot(aes(x = val,
             group = breed)) +
  geom_density(alpha=0.25) +
  coord_flip() +
  facet_wrap(.~factor,
             scales = "free_x",
             ncol = 7,
             nrow = 2) +
  theme_pubr() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


ggsave(filename = paste(workDir,
                        "/gwa/phe/",
                        "cbarq.sample-dists.pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 3 *3,
       height = 2 * 3)

cb.summ %>%
  ggplot() +
  stat_function(aes(color = breed),
                fun = ~ dnorm(.x, score.mean, score.sd))

cb.samp.list = list()
i=1
for (cb in cb.list){
  print(cb)

  cb.phe = anc %>%
    group_by(dog) %>%
    slice_max(pct, n=1) %>%
    mutate(purebred = pct > 0.85) %>%
    mutate(breed = if_else(purebred & pop %in% cbarq$breed,
                           pop,
                           NA_character_)) %>%
    mutate(score = if_else(is.na(breed),
                           rnorm(mean=(cb.all %>% filter(factor == cb) %>% pull(score.mean)),
                                 sd=(cb.all %>% filter(factor == cb) %>% pull(score.sd)),
                                 n=1),
                           rnorm(mean=(cb.summ %>% filter(breed==breed & factor == cb) %>% pull(score.mean)),
                                 sd=(cb.summ %>% filter(breed==breed & factor == cb) %>% pull(score.sd)),
                                 n=1))) %>% 
    ungroup() %>%
    mutate(score = round(score, digits = 2))

  cb.samp.list[[i]] = cb.phe %>%
    mutate(factor=paste("cbarq.breed-agg.",
                        cb,
                        sep = "")) %>%
    select(dog,factor,score)

  write_tsv(x = (cb.phe %>%
                   mutate(fid=dog,
                          iid=dog,
                          phe=score) %>%
                   select(fid,iid,phe)),
            paste(workDir,
                  "/gwa/phe/",
                  "cbarq.breed-agg.",
                  cb,
                  ".tsv",
                  sep = ""),
            col_names = F)

  i=i+1
}

cb.all = rbindlist(cb.samp.list) %>%
  pivot_wider(id_cols = "dog",
              names_from = "factor",
              values_from = "score")

## FOOD ALLERGIES ----

### Perform PCA on tokens ----
docs = Corpus(VectorSource(answers %>%
                             filter(question %in% c(183)) %>%
                             pull(answer)))
docs = tm_map(docs, function(x) stri_replace_all_regex(x, "<.+?>", " "))
docs = tm_map(docs, function(x) stri_replace_all_fixed(x, "t", " "))
docs = tm_map(docs, PlainTextDocument)
docs = tm_map(docs, stripWhitespace)
docs = tm_map(docs, removeWords, stopwords("english"))
docs = tm_map(docs, removePunctuation)
docs = tm_map(docs, tolower)
dtm <- DocumentTermMatrix(docs)
dtm <- as.matrix(dtm)
dim(dtm)
frequency <- colSums(dtm)
frequency <- sort(frequency, decreasing=TRUE)
mots=frequency[frequency>20]
s=dtm[,which(colnames(dtm) %in% names(mots))]
allergen_pca = PCA(s)

allergen_pca_dog =as.data.table(allergen_pca$ind$coord)
allergen_pca_dog$dog = answers %>%
  filter(question %in% c(183)) %>%
  pull(dog)

af.all = allergen_pca_dog %>%
  select(dog,
         Dim.1,
         Dim.2,
         Dim.3,
         Dim.4,
         Dim.5) %>%
  mutate(fid = dog,
         iid = dog)

for (i in seq(1,5)){
  af =paste("Dim.",i,sep="")
  write_tsv((af.all %>% 
               select(fid,iid,phe=!!af) %>% 
               as.data.table()),
            paste(workDir,
                  "/gwa/phe/",
                  "af.",
                  i,
                  ".tsv",
                  sep = ""), 
            col_names = F)
}


### Specific food allergens ----
allergens = as.data.table(s)

allergens$dog = (answers %>%
                   filter(question %in% c(183)) %>%
                   pull(dog))

allergens = allergens %>%
  pivot_longer(cols = -dog,
               names_to = "token",
               values_to = "present")

allergens = as.data.table(allergens)

keywords = c("none" = "none",
             "poultry" = "chicken",
             "poultry" = "turkey",
             "poultry" = "duck",
             "poultry" = "poultry")

for (keyword in keywords){
  
}
allergens[, allergens := 1]
allergens[answer.lower %like% "none", allergens := 0]
allergens[answer.lower == "", allergens := 0]

allergens[, poultry := 0]
allergens[allergens != 0 & answer.lower %like% "chicken", poultry := 1]
allergens[allergens != 0 & answer.lower %like% "turkey", poultry := 1]
allergens[allergens != 0 & answer.lower %like% "duck", poultry := 1]
allergens[allergens != 0 & answer.lower %like% "poultry", poultry := 1]

allergens[, seafood := 0]
allergens[allergens != 0 & answer.lower %like% "fish", seafood := 1]
allergens[allergens != 0 & answer.lower %like% "shellfish", seafood := 1]
allergens[allergens != 0 & answer.lower %like% "salmon", seafood := 1]

allergens[, redmeat := 0]
allergens[allergens != 0 & answer.lower %like% "beef", redmeat := 1]
allergens[allergens != 0 & answer.lower %like% "steak", redmeat := 1]
allergens[allergens != 0 & answer.lower %like% "lamb", redmeat := 1]
allergens[allergens != 0 & answer.lower %like% "buffalo", redmeat := 1]
allergens[allergens != 0 & answer.lower %like% "pork", redmeat := 1]

allergens[, egg := 0]
allergens[allergens != 0 & answer.lower %like% "egg", egg := 1]

allergens[, dairy := 0]
allergens[allergens != 0 & answer.lower %like% "milk", dairy := 1]
allergens[allergens != 0 & answer.lower %like% "dairy", dairy := 1]
allergens[allergens != 0 & answer.lower %like% "yogurt", dairy := 1]
allergens[allergens != 0 & answer.lower %like% "cheese", dairy := 1]

allergens[, grain := 0]
allergens[allergens != 0 & answer.lower %like% "wheat", grain := 1]
allergens[allergens != 0 & answer.lower %like% "corn", grain := 1]
allergens[allergens != 0 & answer.lower %like% "grain", grain := 1]

allergens[, any := 0]
allergens[poultry == 1 | seafood == 1 | redmeat == 1 | egg == 1 | dairy == 1 | grain == 1, any := 1]

## BREED ANCESTRY for ----
anc %>%
  group_by(pop) %>%
  summarize(sum_pct=sum(pct),
            mean_pct=mean(pct),
            sd_pct=sd(pct)) %>%
  ungroup() %>%
  mutate(tot_prop = sum_pct/sum(sum_pct)) %>%
  arrange(-tot_prop) %>%
  head(n=5)

## Breed aesthetics CCA (do last because library for CCA fucks up dplyr) ----
# library(CCA)
# mq.cca.input = answers %>%
#   filter(question %in% c(121,122,123,124,125,126,127,128)) %>%
#   mutate(trait = paste("Q", question, sep = "")) %>%
#   mutate(answer = as.numeric(answer)) %>%
#   dplyr::select(dog,trait,answer) %>%
#   pivot_wider(id_cols = "dog",
#               names_from = "trait",
#               values_from = "answer", 
#               values_fn = min) %>%
#   merge((dogs %>% dplyr::select(dog,breed)))
# 
# mq.cca.dogs = mq.cca.input$dog
# 
# mq.cca.breeds = data.table(breed=dogs$breed) %>%
#   group_by(breed) %>%
#   summarize(n=n()) %>%
#   filter(!is.na(breed)) %>%
#   filter(n>100) %>%
#   pull(breed)
# 
# y = mq.cca.input %>%
#   filter(!is.na(breed)) %>%
#   mutate(count=1) %>%
#   pivot_wider(names_from = "breed",
#               values_from = "count",
#               values_fn = min, 
#               values_fill = 0) %>%
#   dplyr::select(all_of(c("dog",mq.cca.breeds))) %>%
#   group_by(dog) %>%
#   slice_head(n=1) %>%
#   ungroup() %>%
#   arrange(dog)
# 
# x = mq.cca.input %>%
#   filter(dog %in% y$dog) %>%
#   group_by(dog) %>%
#   slice_head(n=1) %>%
#   ungroup() %>%
#   arrange(dog)
# 
# y = y %>%
#   filter(dog %in% x$dog) %>%
#   arrange(dog)
# 
# x = x %>%
#   filter(dog %in% y$dog) %>%
#   arrange(dog)
# 
# y.breeds = colnames(y)[-1]
# x.traits = colnames(x)[-1]
# 
# x.in = as.matrix((x %>% 
#                     dplyr::select(-dog) %>% 
#                     dplyr::select(-breed)))
# y.in = as.matrix((y %>% 
#                     dplyr::select(-dog)))
# 
# dim(x.in)
# dim(y.in)
# 
# mq.cca = cc(X = x.in,
#             Y = y.in)
# 
# mq.cca.scores = x.in %*% t(as.matrix(mq.cca$ycoef))
# mq.cca.scores = as.data.frame(mq.cca.scores)
# mq.cca.scores$dog = x$dog
# 
# mq.cca.scores %>%
#   merge((dogs %>% dplyr::select(dog,breed)), by = "dog") %>%
#   filter(breed %in% mq.cca.breeds) %>%
#   group_by(breed) %>%
#   summarise(across(everything(), mean, na.rm=TRUE)) %>%
#   View()

# SAVE DATA ----
## Finalize tables ----


## Save files ----
write_tsv(bq.all,
          paste(
            workDir,
            "/gwa/",
            "DarwinsArk_",
            "bq.all",
            ".tsv",
            sep = ""
          ))

write_tsv(mq.all,
          paste(
            workDir,
            "/gwa/",
            "DarwinsArk_",
            "mq.all",
            ".tsv",
            sep = ""
          ))

write_tsv(mp.all,
          paste(
            workDir,
            "/gwa/",
            "DarwinsArk_",
            "mp.all",
            ".tsv",
            sep = ""
          ))

write_tsv(fa.agg.all,
          paste(
            workDir,
            "/gwa/",
            "DarwinsArk_",
            "fa.breed-agg.all",
            ".tsv",
            sep = ""
          ))

write_tsv(cb.all,
          paste(
            workDir,
            "/gwa/",
            "DarwinsArk_",
            "cbarq.breed-agg.all",
            ".tsv",
            sep = ""
          ))

## Save image ----
save.image(file = paste(workDir,
                        "/dat/",
                        paste(
                          paste("DarwinsArk",
                                freezeDate,
                                "phenotyping",
                                sep = "_"),
                          ".RData",
                          sep = ""
                        ),
                        sep = ""))

## Define combos ----
all.pheno = dogs %>%
  select(dog,sex,age,
         household_cats,household_dogs,
         drug_psych,drug_allergy,drug_pain,drug_sedative) %>%
  merge(mq.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(mp.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(bq.all %>% unique(),
        by.x = "dog",
        by.y = "iid",
        all = T) %>%
  merge(fa.all %>% unique(),
        by.x = "dog",
        by.y = "iid",
        all = T) %>%
  merge(lineage.ext %>% unique(),
        by = "dog",
        all = T) %>%
  merge(lineage.int %>% unique(),
        by = "dog",
        all = T) %>%
  merge(fa.agg.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(cb.all %>% unique(),
        by = "dog",
        all = T) %>%
  group_by(dog) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  dplyr::select(-c("fid.x","fid.y")) %>%
  filter(dog %in% gwas$dog)

write_tsv((all.pheno %>%
             mutate(FID=dog,
                    IID=dog,
                    .before = dog) %>%
             select(-dog)),
          paste(workDir,
                "/gwa/phe/",
                "all.tsv", sep = ""), col_names = F)

write_lines((all.pheno %>%
               select(-dog) %>%
               colnames()),
            paste(workDir,
                  "/gwa/phe/",
                  "/all.key.txt", 
                  sep = ""))

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

all.pairs = expand.grid.unique(x = (all.pheno %>%
                                      select(-dog) %>%
                                      colnames()),
                               y = (all.pheno %>%
                                      select(-dog) %>%
                                      colnames())) %>%
  as.data.frame() %>%
  rename(A = V1,
         B = V2) %>%
  rowwise() %>%
  mutate(pA = which(colnames(all.pheno)[-1] == A),
         pB = which(colnames(all.pheno)[-1] == B))

write_tsv(all.pairs,
          paste(workDir,
                "/gwa/phe/",
                "/all.pairs.tsv", 
                sep = ""), col_names = F)

# Specific splits:

cb.pairs = all.pairs %>%
  filter(B %like% "cbarq") %>%
  filter((A %like% "lineage") | (A %like% "-mean-filled") | (A %like% "mq.") | (A %like% "mp."))

write_tsv(cb.pairs,
          paste(workDir,
                "/gwa/phe/",
                "/cb.pairs.tsv", 
                sep = ""), col_names = F)

#### Focused pairs: ----
cbarq.vs.all = cba
all.pheno = dogs %>%
  select(dog,sex,neutered,purebred,age_mean,
         environ_filled,origin_filled,
         household_scale_ppl,household_cats,household_dogs,
         train_level,drug_psych,drug_allergy,drug_pain,drug_sedative) %>%
  merge(mq.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(mp.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(af.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(bq.all %>% unique(),
        by = "dog",
        all = T) %>%
  merge(fa.all %>% unique(),
        by = "dog",
        all = T) %>%
  group_by(dog) %>%
  slice_head(n=1) %>%
  ungroup()

write_tsv((all.pheno %>%
             mutate(FID=dog,
                    IID=dog,
                    .before = dog) %>%
             select(-dog)),
          paste("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/data/pheno/",freezeDate,"/all.tsv", sep = ""), col_names = F)

write_lines((all.pheno %>%
               select(-dog) %>%
               colnames()),
            paste("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/data/pheno/",freezeDate,"/all.key.txt", sep = ""))

expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

all.pairs = expand.grid.unique(x = (all.pheno %>%
                                      select(-dog) %>%
                                      colnames()),
                               y = (all.pheno %>%
                                      select(-dog) %>%
                                      colnames())) %>%
  as.data.frame() %>%
  rename(A = V1,
         B = V2) %>%
  rowwise() %>%
  mutate(pA = which(colnames(all.pheno)[-1] == A),
         pB = which(colnames(all.pheno)[-1] == B))

write_tsv(all.pairs,
          paste("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/data/pheno/",freezeDate,"/all.pairs.tsv", sep = ""), col_names = F)

# MORE PLOTS ----

## Factor Scores ----
library(ggstatsplot)

for (fa in unique(scores$factor)){

fa.df = scores %>%
  filter(factor == fa) %>%
  merge(dogs, by = "dog", all.x = T) %>%
  mutate(gwas = dog %in% gwas$dog) %>%
  group_by(gwas) %>%
  mutate(n = n()) %>%
  mutate(set = if_else(gwas,
                       paste("genetic data\n","(n=",format(n,nsmall=1, big.mark=","),")", sep =""),
                       paste("survey data\n","(n=",format(n,nsmall=1, big.mark=","),")", sep =""))) %>%
  ungroup()

quants = fa.df %>% 
  pull(score_fill_samp) %>% 
  quantile()

quants.pal = c("#028AA5",
               "#3D7086",
               "#775566",
               "#B23B47",
               "#EC2027")

fa.df = fa.df %>%
  mutate(quartile = ntile(score_fill_samp, 5)) %>%
  mutate(color = quants.pal[quartile])

p = fa.df %>%
  ggplot(aes(x = set,
             y = score_fill_samp)) +
  coord_flip() +
  geom_jitter(aes(color = color),
             alpha = 0.25,width = 0.25) +
  geom_violin(fill = "#FFFFFF",
              alpha = 0.50) +
  geom_boxplot(fill = "#FFFFFF",
               color = "#505050",
               alpha = 0.50, 
               notch = T,
               width = 0.5,
               outlier.shape = NA) +
  scale_color_identity() +
  theme_pubr() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())


ggsave(filename = paste(workDir,
                        "/gwa/phe/",
                        "density-plot.",
                        fa,
                        ".pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       dpi = 300,
       units = "in",
       width = 3*2,
       height = 1.5*2)
}


### PLOT ANY ----
q=36
bq="bq.036"

thisQuestion = answers %>%
  filter(question == q) %>%
  select(question,dog,answer,option,age) %>%
  mutate(answer = if_else(option %in% c("I don't know",
                                        "I'm not sure",
                                        "Maybe"),
                          "",
                          answer)) %>%
  mutate(fid = dog,
         iid = dog,
         phe = as.numeric(answer)) %>%
  merge((dogs %>%
           select(dog,
                  hgt=height_filled)),
        by = "dog", all.x = T)



p = thisQuestion %>%
                 filter(!is.na(phe)) %>%
                 filter(dog %in% gwas$dog) %>%
                 mutate(n = n()) %>%
                 mutate(set = paste("genetic data ",
                                    "(n= ",n,")",
                                    sep = "")) %>%
                 select(dog,phe,set,option) %>%
  mutate(option = if_else(option == "Neither Agree Nor Disagree",
                          "Neither",
                          option)) %>%
  mutate(lab = paste(as.character(as.integer(phe)),
                     ". ", option,
                     sep = "")) %>%
  ungroup() %>%
  mutate(mean_phe = mean(as.integer(phe), na.rm = T)) %>%
  group_by(phe,lab,set,mean_phe) %>%
  summarize(phe_n = n()) %>%
  mutate(mean_opt = if_else(phe == round(mean_phe),
                            lab,
                            NA_character_)) %>%
  ggplot(aes(x = lab,
             y = phe_n)) +
  coord_flip(expand = F) +
  geom_col(linewidth = 0.75, width = 0.5) +
  labs(x = "option",
       y = "# dogs") +
  theme_pubr() +
  theme(legend.position="bottom", 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=45,
                                   hjust = 1),
        axis.text.y = element_text(size=20),
        plot.margin = margin(0,0.25,0,0, "cm"),
        plot.title = element_text(size=10))

ggsave(filename = paste(workDir,
                        "/gwa/phe/",
                        "dist-plot.ppt.",
                        bq,
                        ".png",
                        sep = ""),
       plot = p,
       device = "png",
       dpi = 300,
       units = "in",
       width = 3*1.25,
       height = 3.5*1.25)

