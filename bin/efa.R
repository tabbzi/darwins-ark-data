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
require(wesanderson)
require(jtools)

### Correlation and Regression ----
require(pspearman)
require(jtools)
require(caret)

### Dimensional Reduction ----
require(psych)
require(nFactors)
require(FactoMineR)
require(factoextra)
require(corrr)
require(car)
require(BBmisc)
require(mice)

## ARGUMENTS ----
parser <- ArgumentParser()
parser$add_argument("--dir",
                    help = "Working directory", 
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data")
parser$add_argument("--date",
                    help = "Data freeze date",
                    default = "20221120")

args <- parser$parse_args()
print(args)

### Set working directory ----
workDir = args$dir
setwd(workDir)

### Set data freeze ----
freezeDate = args$date
freezeDir = paste(workDir,"data","database","data_freeze",freezeDate, sep = "/")

## LOAD DATA ----
dogs = read_csv(paste(workDir,
                      "/dat/",
                      "DarwinsArk_",
                      as.character(freezeDate),
                      "_dogs.csv",
                      sep=""))

questions = read_csv(paste(workDir,
                           "/dat/",
                           "DarwinsArk_",
                           as.character(freezeDate),
                           "_questions.csv",
                           sep=""))

answers = read_csv(paste(workDir,
                         "/dat/",
                         "DarwinsArk_",
                         as.character(freezeDate),
                         "_answers.csv",
                         sep=""))

# EXPLORATORY FACTOR ANALYSIS ----

## Select ordinal survey items ----
q_fa = questions %>%
  filter(format %in% c(1,2,17,19,20,25) & id <= 250) %>%
  pull(id)

## Pull responses ----
data = answers %>%
  filter(question %in% q_fa) %>%
  filter(!is.na(answer)) %>%
  select(question,dog,answer) 

## Scale responses and generate matrix ----
raw_mat = data %>%
  group_by(question) %>%
  mutate(ans_norm = normalize(as.numeric(answer), method = "standardize", range = c(0,6))) %>%
  ungroup() %>%
  pivot_wider(id_cols = "dog",
              names_from = "question",
              values_from = "ans_norm")

## Drop dogs with any missing data ----
mat_na_drop = raw_mat %>% drop_na()
mat_dogs = mat_na_drop$dog
mat = mat_na_drop %>% select(-dog)
mat_ques = colnames(mat)

## Get counts ----
qN = ncol(mat)
dN = nrow(mat)

## Set output ----
output = paste(workDir,
               "/dat/DarwinsArk_",
               freezeDate,
               "_factor-analysis",
               "_discovery-no-na",
               "_qn-",
               qN,
               "_dn-",
               dN,
               sep = "")

write_lines(x = paste(qN,"items"),
            file = paste(output,
                         ".metadata.txt",
                         sep = ""),
            append = T)

write_lines(x = paste(dN,"dogs"),
            file = paste(output,
                         ".metadata.txt",
                         sep = ""),
            append = T)

## Generate correlation matrix ----
mat_cor <- cor(as.matrix(mat))

## Kaiser-Meyer-Olkin (KMO) test ----
fa.kmo = KMO(r=mat_cor) # Overall MSA = 0.94

write_lines(x = paste("KMO Overall MSA:", fa.kmo$MSA),
            file = paste(output,
                         ".metadata.txt",
                         sep = ""),
            append = T)

## Bartlett’s Test of Sphericity ----
fa.bart = cortest.bartlett(mat_cor, n = nrow(mat)) # p = 0
write_lines(x = paste("Bartlett Test of Sphericity: p=", fa.bart$p),
            file = paste(output,
                         ".metadata.txt",
                         sep = ""),
            append = T)

## Check whether determinant is + or - ----
fa.det = det(cor(mat_cor))

write_lines(x = paste("Determinant =", fa.det),
            file = paste(output,
                         ".metadata.txt",
                         sep = ""),
            append = T)

## Perform parallel analysis ----
fa.ap <- parallel(subject=nrow(mat),
                  var=ncol(mat),
                  rep=100,
                  cent=.05)

## Perform scree test ----
fa.ev <- eigen(cor(as.matrix(mat))) 
fa.ns <- nScree(x=fa.ev$values,
                aparallel=fa.ap$eigen$qevpea)

## Perform nfactors test ----
fa.nf = nfactors(x = mat,
                        n = 50,
                        fm = "pa",
                        rotate = "varimax")

## Plot heuristic test results ----
pdf(file = paste(output,
                     ".scree-test.pdf",
                     sep = ""),
    width = 6,
    height = 6)

plotnScree(fa.ns)

dev.off()

pdf(file = paste(output,
                 ".n-factors.pdf",
                 sep = ""),
    width = 6,
    height = 6)

nfactors(x = mat,
         n = 50,
         fm = "pa",
         rotate = "varimax")

dev.off()

## Summarize heuristic results ----
bind_cols(fa.ns$Analysis,
          fa.nf$vss.stats)

## Set number of factors to extract ----
nFac = 25
write_lines(x = paste("Factors to extract =", nFac),
            file = paste(output,
                         ".metadata.txt",
                         sep = ""),
            append = T)

## Perform exploratory factor analysis ----
fa.varimax <- fa(r=mat_cor,
                 nfactors = nFac,
                 fm="pa",
                 max.iter=100,
                 rotate="varimax")

# Item Uniqueness: uniqueness, sometimes referred to as noise, corresponds to the proportion of variability, which can not be explained by a linear combination of the factors. A high uniqueness for a variable indicates that the factors do not account well for its variance. (what exists in the unique variance instead of common variance)
fa.varimax.uniqueness = data.table(question = names(fa.varimax$uniquenesses),
                                   uniqueness = fa.varimax$uniquenesses) %>%
  merge((questions %>% arrange(id) %>% select(id,string,source=tags,options,survey=title) %>% mutate(id = as.character(id))), by.x = "question", by.y = "id") %>% arrange(as.numeric(question))

write_tsv(fa.varimax.uniqueness,
          paste(output,
                "_fn-",
                nFac,
                ".uniqueness.tsv",
                sep = ""))

# Factor Loadings
fa.varimax.str = fa.varimax$loadings
fa.varimax.pat = promax(fa.varimax$loadings)$loadings

fa.varimax.loadings = list()
for (i in 1:nFac){
  load_pat = fa.varimax.pat[abs(fa.varimax.str[,i]) > 0.2,i]
  load_str = fa.varimax.str[abs(fa.varimax.str[,i]) > 0.2,i]
  fa.varimax.loadings[[i]] = merge(data.table(factor = colnames(fa.varimax.str)[i],
                                              question = as.numeric(names(load_pat)),
                                              pattern = load_pat),
                                   data.table(factor = colnames(fa.varimax.str)[i],
                                              question = as.numeric(names(load_str)),
                                              structure = load_str),
                                   by = c("factor","question"),
                                   all = T)
}

fa.varimax.loadings = rbindlist(fa.varimax.loadings) %>%
  merge((questions %>% select(question=id,text=string,options,source=tags)), by = "question", all.x = T) %>%
  select(factor,question,text,source,options,pattern,structure) %>%
  mutate(factor = as.numeric(gsub("PA","",factor))) %>%
  arrange(factor,-abs(pattern))

fa.varimax.loadings = fa.varimax.loadings %>%
  mutate(options = if_else(options == "Strongly Agree|Agree|Neither Agree Nor Disagree|Disagree|Strongly Disagree",
                           "Strongly Agree to Strongly Disagree",
                           if_else(options == "Never|Rarely|Sometimes|Often|Always",
                                   "Never to Always",
                                   if_else(options == "Not at all true|Somewhat true|Mainly true|Definitely true",
                                           "Not True to Definitely True",
                                           if_else(options == "None|Less than 30 min|30 min - 1 hr|1 - 2 hr|2 - 3 hr|3 - 5 hr|5+ hr",
                                                   "None to 5+ hr",
                                                   if_else(options == "Rarely|Once a month|Once a week|Once a day or more",
                                                           "Rarely to Once a Day",
                                                           options))))))

write_tsv(fa.varimax.loadings,
          paste(output,
                "_fn-",
                nFac,
                ".loadings.tsv",
                sep = ""))

fa.varimax.scores = factor.scores(x = mat, f = fa.varimax)
fa.varimax.scores.parsed = fa.varimax.scores$scores %>% as.data.table()
fa.varimax.scores.parsed$dog = mat_dogs
fa.varimax.scores.parsed = fa.varimax.scores.parsed %>%
  select(dog,paste("PA", 1:nFac, sep = ""))

# Summary Statistics
fa.question.stats = data.table(describe(mat)) %>%
  mutate(question = mat_ques)

fa.question.stats = fa.question.stats %>%
  merge((questions %>% mutate(id = as.character(id))), by.x = "question", by.y = "id", all.x = T, all.y = T) %>%
  mutate(included = question %in% fa.varimax.loadings$question) %>%
  arrange(as.numeric(question))

write_tsv(fa.question.stats,
          paste(output,
                "_fn-",
                nFac,
                ".stats-items.tsv",
                sep = ""))

Lambda <- unclass(fa.varimax$loadings)
p <- nrow(Lambda)
f <- ncol(Lambda)

vx <- colSums((fa.varimax$loadings)^2)
varex <- rbind(`SS loadings` = vx)

if (is.null(attr(fa.varimax$loadings, "covariance"))) {
  varex <- rbind(varex, `Proportion Var` = vx/p)
  if (f > 1) 
    varex <- rbind(varex, `Cumulative Var` = cumsum(vx/p))
}

fa.factor.stats = tibble::rownames_to_column(as.data.frame(varex), "x")

fa.factor.stats = fa.factor.stats %>%
  pivot_longer(cols = paste("PA", seq(1,nFac), sep = ""),
               names_to = "factor",
               names_prefix = "PA",
               values_to = "statistic") %>%
  pivot_wider(names_from = "x",
              values_from = "statistic") %>%
  select(Factor = factor, `Sum of Square Loadings` = `SS loadings`, `Proportion Variance` = `Proportion Var`, `Cumulative Variance` = `Cumulative Var`) %>%
  arrange(as.numeric(Factor))

write_tsv(fa.factor.stats,
          paste(output,
                "_fn-",
                nFac,
                ".stats-factors.tsv",
                sep = ""))

require(corrplot)
cx = cor(fa.varimax$loadings)
rownames(cx)= gsub("PA","",rownames(cx))
colnames(cx) = rownames(cx)
cx <- cx[order(as.integer(rownames(cx))), 
         order(as.integer(colnames(cx)))]
pdf(paste(output,
          "_fn-",
          nFac,
          ".corr-plot.pdf",
          sep = ""))
corrplot.mixed(cx, lower.col = "grey", number.cex = 0.4, mar = c(0,0,0,2))
dev.off()

rowcx = rownames(cx)
colcx = colnames(cx)

fa.cors = cx %>% 
  as_tibble() %>% 
  mutate(fac1 = rowcx) %>% 
  pivot_longer(cols = c(colcx,-fac1), 
               values_to = "cor", 
               names_to = "fac2") %>% 
  filter(fac1!=fac2) %>% 
  arrange(-abs(cor))

write_tsv(fa.cors,
          paste(output,
                "_fn-",
                nFac,
                ".factor-corr.tsv",
                sep = ""))

# ITEM CORRELATIONS ----
item.cx = cor(mat)
rownames(item.cx)= mat_ques
colnames(item.cx) = mat_ques
item.cx <- item.cx[order(as.integer(rownames(item.cx))), order(as.integer(colnames(item.cx)))]

## Plot item-item correlations ----
pdf(file = paste(output,
                 ".items",
                 ".corr-plot.pdf",
                 sep = ""),
    width = 24,
    height = 24)
corrplot.mixed(item.cx,
               lower.col = "grey",
               number.cex = 0.4,
               mar = c(0,0,0,2))
dev.off()

## Plot intra-survey item-item correlations ----
for (S in (questions %>% 
           filter(id %in% mat_ques) %>% 
           pull(survey) %>% 
           unique())){
  p = mat %>%
    corrr::correlate() %>%
    corrr::focus((questions %>% filter(survey == S) %>% filter(id %in% mat_ques) %>% arrange(as.integer(id)) %>% pull(id) %>% as.character()), mirror = TRUE) %>% 
    corrr::network_plot(colors = c("red", "green"))
  ggsave(filename = paste(output,
                          ".items",
                          ".corr.",
                          "survey-",
                          S,
                          ".pdf",
                          sep = ""),
         plot = p,
         device = "pdf",
         width = 8,
         height = 8,
         units = "in")
}

## Generate table of item-item correlations ----
item.cx.long = item.cx %>%
  as_cordf() %>%
  stretch() %>%
  filter(!is.na(r) & x!=y) %>%
  merge((questions %>% select(x_item = id, x_string = string, x_options = options)), by = "x", by.y = "x_item", all.x = T) %>%
  merge((questions %>% select(y_item = id, y_string = string, y_options = options)), by = "y", by.y = "y_item", all.x = T)

## Save table ----
write_tsv(x = item.cx.long,
          file = paste(output,
                       ".items",
                       ".corr.tsv",
                       sep = ""))

# FACTOR SCORING ----
## Set seed for consistent sampling ----
set.seed(seed = 101)

## Function to replace NAs with random sample ----
replace_func <- function(x, y) {
  inds <- is.na(x)
  if (length(x) > 1 & any(inds)) {
    x[inds] <- sample(x[!inds], sum(inds), replace = T)
    x
  }
  else if(any(inds)) {
    x[inds] <- sample(y[!is.na(y)], 1, replace = T)
    x
  } else x
}

## Fill missing data and calculate scores ----
fa.factor.scores.filled = list()

for (thisFac in fa.factor.stats$Factor){
  qfac = fa.varimax.loadings %>%
    filter(factor == thisFac) %>%
    pull(question)
  
  qfacscore = answers %>%
    filter(question %in% qfac) %>%
    filter(!is.na(answer)) %>%
    select(question,
           dog,
           ans = answer,
           age) %>%
    as_tibble() %>%
    mutate(ans = as.numeric(ans)) %>%
    mutate(original = TRUE) %>%
    complete(dog,question) %>%
    group_by(question) %>%
    mutate(ans_mean = mean(ans, na.rm = T),
           ans_median = median(ans, na.rm = T)) %>%
    ungroup() %>%
    group_by(dog) %>%
    mutate(total = n(),
           nomiss = sum(original, na.rm = T)) %>%
    mutate(completeness = nomiss/total) %>%
    group_by(question) %>%
    mutate(ans_real = ans) %>%
    mutate(ans_fill_samp = replace_func(ans, .$ans),
           ans_fill_mean = if_else(is.na(ans),
                                   ans_mean,
                                   ans)) %>%
    group_by(question) %>%
    mutate(ans_real_norm = normalize(ans_real, 
                                     method = "standardize", 
                                     range = c(0,6)),
           ans_fill_samp_norm = normalize(ans_fill_samp, 
                                     method = "standardize", 
                                     range = c(0,6)),
           ans_fill_mean_norm = normalize( ans_fill_mean, 
                                     method = "standardize", 
                                     range = c(0,6))) %>%
    ungroup() %>%
    merge((fa.varimax.loadings %>%
             filter(factor == thisFac) %>%
             select(question,structure)),
          by = "question",
          all.x = T) %>%
    mutate(product_real = structure*ans_real_norm,
            product_fill_samp = structure*ans_fill_samp_norm,
           product_fill_mean = structure*ans_fill_mean_norm) %>%
    group_by(dog) %>%
    summarize(score_real = sum(product_real),
              score_fill_samp = sum(product_fill_samp),
              score_fill_mean = sum(product_fill_mean),
              completeness = mean(completeness),
              age_mean = mean(age, na.rm = T)) %>%
    filter(completeness >= 0.75) %>%
    ungroup() %>%
    mutate(score_real_norm = scale(score_real),
           score_fill_samp_norm = scale(score_fill_samp),
           score_fill_mean_norm = scale(score_fill_mean))
  
  qfacscore$factor = thisFac
  fa.factor.scores.filled[[thisFac]] = qfacscore
}

fa.factor.scores.filled = rbindlist(fa.factor.scores.filled)

# SAVE DATA ----
## Finalize tables ----
fa.factor.scores = fa.factor.scores.filled %>%
  merge((fa.varimax.scores.parsed %>%
           pivot_longer(cols = starts_with("PA"),
                        names_to = "factor",
                        names_prefix = "PA",
                        values_to = "score_raw") %>%
           mutate(discovery = TRUE)),
        by = c("dog","factor"),
        all.x = T) %>%
  mutate(factor = as.integer(factor)) %>%
  arrange(discovery,dog,factor) %>%
  mutate(discovery = if_else(discovery == TRUE,
                             TRUE,
                             FALSE,
                             FALSE)) %>%
  select(factor,
         dog,
         discovery,
         age_mean,
         score_raw,
         score_real,
         score_real_norm,
         score_fill_samp,
         score_fill_samp_norm,
         score_fill_mean,
         score_fill_mean_norm)

fa.factor.scores %>% 
  ggplot(aes(x = score_raw, 
             y = score_real_norm)) + 
  geom_point(alpha=0.2) + 
  theme_pubr()

## Save files ----
write_csv(fa.factor.scores,
          paste(output,
                ".factor-scores.csv",
                sep=""))

## Save image ----
save.image(file = paste(workDir,
                        "/dat/", 
                        paste(paste("DarwinsArk",
                                    freezeDate,
                                    "factor-analysis",
                                    sep="_"), 
                              ".RData", 
                              sep=""), 
                        sep=""))
