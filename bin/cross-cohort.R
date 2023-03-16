library(tidyverse)
library(data.table)
library(haven)
library(sjlabelled)
library(measurements)

# Load Darwin's Ark data: ----
darwinsark.dogs = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_dogs.csv")

darwinsark.questions = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_questions.csv")

darwinsark.answers = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_answers.csv")

darwinsark.geno = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/dat/gen/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465.fam", col_names = c("fid","dog","mat","pat","sex","phe"))

# Load Dog Aging Project data: ----
load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_HLES_dog_owner_v1.0.RData")
load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_AFUS_dog_owner_v1.0.RData")
load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_CSLB_v1.0.RData")

# Define common survey items ----
trait.map = data.table(phecode.darwinsark = 
                         c("bq.206",
                           "bq.227",
                           "bq.207",
                           "bq.208",
                           "bq.209",
                           "bq.210",
                           "bq.211",
                           "bq.240",
                           "bq.212",
                           "bq.228",
                           "bq.213",
                           "bq.214",
                           "bq.215",
                           "bq.216",
                           "bq.229",
                           "bq.217",
                           "bq.230",
                           "bq.231",
                           "bq.218",
                           "bq.219",
                           "bq.220",
                           "bq.221",
                           "bq.232",
                           "bq.222",
                           "bq.223",
                           "bq.233",
                           "bq.234",
                           "bq.235",
                           "bq.224",
                           "bq.225",
                           "bq.226",
                           "bq.237",
                           "bq.238",
                           "bq.239",
                           "bq.037",
                           "bq.061",
                           "bq.093",
                           "bq.080",
                           "bq.095",
                           "bq.036"),
                       phecode.dap = 
                         c("afus_dora_1_excited_food",
                           "afus_dora_2_most_walks_off_leash",
                           "afus_dora_3_human_leftovers_in_bowl",
                           "afus_dora_4_waits_for_scraps",
                           "afus_dora_5_choosy_treats",
                           "afus_dora_6_waits_during_food_prep",
                           "afus_dora_7_turns_down_food",
                           "afus_dora_8_finishes_meal_quickly",
                           "afus_dora_9_inspects_unfamilar_foods",
                           "afus_dora_10_runs_around_alot",
                           "afus_dora_11_interested_eating_after_meal",
                           "afus_dora_12_slow_eater",
                           "afus_dora_13_eats_treats_quickly",
                           "afus_dora_14_human_food_during_meals",
                           "afus_dora_15_eat_anything",
                           "afus_dora_16_very_fit",
                           "afus_dora_17_human_food_often",
                           "afus_dora_18_upset_stomach_some_foods",
                           "afus_dora_19_should_lose_weight",
                           "afus_dora_20_most_walks_on_leash",
                           "afus_dora_21_restrict_exercise",
                           "afus_dora_22_diet_to_control_weight",
                           "afus_dora_23_hungry_all_time",
                           "afus_dora_24_walks_high_energy",
                           "afus_dora_25_careful_dogs_weight",
                           "afus_dora_26_sensitive_stomach",
                           "afus_dora_27_very_greedy",
                           "afus_dora_28_vet_often_for_health_problems",
                           "afus_dora_29_happy_dogs_weight",
                           "afus_dora_31_regulate_exercise_for_weight",
                           "afus_dora_32_gets_lots_of_exercise",
                           "afus_dora_33_upset_stomach_often",
                           "afus_dora_34_no_human_food_during_meals",
                           "afus_dora_35_eat_non_food_objects",
                           "cslb_avoid",
                           "cslb_find_food",
                           "cslb_pace",
                           "cslb_recognize",
                           "cslb_stare",
                           "cslb_stuck"))

# Get DORA and CCDR from Darwin's Ark data ----
darwinsark.common = darwinsark.answers %>%
  filter(question %in% (darwinsark.questions %>%
                          filter(tags %in% c("DORA","CCDR")) %>%
                          pull(id))) %>%
  mutate(item = paste("bq.",str_pad(question,width = 3,side = "left","0"),sep=""),
         value = as.numeric(answer)) %>%
  filter(item %in% trait.map$phecode.darwinsark) %>%
  select(dog,item,value) %>%
  mutate(cohort = "Darwin's Ark") %>%
  merge(trait.map, 
        by.x = "item", 
        by.y = "phecode.darwinsark") %>%
  rename(match = phecode.dap)

## Define levels ----
darwinsark.common.levels = darwinsark.common %>% 
  group_by(item) %>%
  select(item,match,value) %>%
  unique() %>%
  arrange(item,value) %>%
  mutate(id = row_number()) %>%
  ungroup()

# Get DORA and CCDR from Dog Aging Project data ----
dogagingproject.common = bind_rows((AFUS_dog_owner %>%
                                      select(dog = dog_id,colnames(AFUS_dog_owner)[colnames(AFUS_dog_owner) %in% trait.map$phecode.dap]) %>%
                                      pivot_longer(cols = -dog,
                                                   names_to = "item",
                                                   values_to = "value") %>%
                                      mutate(value = as.numeric(value)) %>%
                                      mutate(cohort = "Dog Aging Project") %>%
                                      merge(trait.map, 
                                            by.x = "item", 
                                            by.y = "phecode.dap") %>%
                                      rename(match = phecode.darwinsark)),
                                   (CSLB %>%
                                      filter(cslb_year == 2020) %>%
                                      select(dog = dog_id,colnames(CSLB)[colnames(CSLB) %in% trait.map$phecode.dap]) %>%
                                      pivot_longer(cols = -dog,
                                                   names_to = "item",
                                                   values_to = "value") %>%
                                      mutate(value = as.numeric(value)) %>%
                                      mutate(cohort = "Dog Aging Project") %>%
                                      merge(trait.map, 
                                            by.x = "item", 
                                            by.y = "phecode.dap") %>%
                                      rename(match = phecode.darwinsark)))

dogagingproject.common.levels = dogagingproject.common %>% 
  group_by(item) %>%
  select(item,match,value) %>%
  unique() %>%
  arrange(item,value) %>%
  mutate(id = row_number()) %>%
  ungroup()

# Identify corresponding values ----
item.levels = merge(darwinsark.common.levels,
                    dogagingproject.common.levels,
                    by.x= c("item","match","id"),
                    by.y= c("match","item","id"),
                    all = T) %>%
  na.omit()

# Recode Darwin's Ark data ----
darwinsark.merge = darwinsark.common %>%
  merge(item.levels, 
        by.x = c("item","match","value"),
        by.y = c("item","match","value.x"),
        all = T)

# Merge together ----
all.dogs = bind_rows((darwinsark.merge %>%
                        select(darwinsark.item = item,
                               dogagingproject.item = match,
                               cohort, 
                               dog, 
                               value=value.y)),
                     (dogagingproject.common %>%
                        select(darwinsark.item = match,
                               dogagingproject.item = item,
                               cohort,
                               dog,
                               value)))

# Extract merged phenotypes ----
all.dogs = all.dogs %>%
  mutate(item = paste(darwinsark.item,dogagingproject.item,sep="|")) %>%
  group_by(item) %>%
  mutate(val.norm = scale(value))

all.dogs %>%
  filter(item == "bq.036|cslb_stuck") %>%
  ggplot(aes(x = value,
             fill = cohort)) +
  geom_bar(position = "dodge") +
  labs(title="bq.036|cslb_stuck") +
  theme_pubr()

all.dogs.pheno = all.dogs %>%
  ungroup() %>%
  mutate(item = darwinsark.item,
         fid = if_else(cohort == "Dog Aging Project",
                       "dogagingproject",
                       "darwinsark"),
         iid = dog,
         phe = value) %>%
  arrange(fid,iid) %>%
  select(item,fid,iid,phe)

for (P in unique(all.dogs.pheno$item)){
  write_tsv(all.dogs.pheno %>% filter(item==P) %>% select(-item),
            paste("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/2023-01-01/cross-cohort.DORA.CCDR.",P,".tsv",sep=""),
            col_names = F)
}

# Factor analysis ----
require(psych)
require(nFactors)
require(FactoMineR)
require(factoextra)
require(corrr)
require(car)
require(BBmisc)
require(mice)

## Scale responses and generate matrix ----
raw_mat = all.dogs %>%
  mutate(item.uniq = paste(darwinsark.item,dogagingproject.item,sep="|")) %>%
  mutate(dog.uniq = paste(cohort,dog,sep="|")) %>%
  select(-c("darwinsark.item",
            "dogagingproject.item",
            "dog",
            "cohort")) %>%
  pivot_wider(id_cols = "dog.uniq",
              names_from = "item.uniq",
              values_from = "value",
              values_fn = min)

## Drop dogs with any missing data ----
mat_na_drop = raw_mat %>% drop_na()
mat_dogs = mat_na_drop$dog.uniq
mat = mat_na_drop %>% select(-dog.uniq)
mat_ques = colnames(mat)

## Get counts ----
qN = ncol(mat)
dN = nrow(mat)

## Generate correlation matrix ----
mat_cor <- cor(as.matrix(mat))

## Kaiser-Meyer-Olkin (KMO) test ----
fa.kmo = KMO(r=mat_cor) # Overall MSA = 0.82
fa.kmo

## Bartlett’s Test of Sphericity ----
fa.bart = cortest.bartlett(mat_cor, n = nrow(mat)) # p = 0
fa.bart

## Check whether determinant is + or - ----
fa.det = det(cor(mat_cor))
fa.det

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
plotnScree(fa.ns)

nfactors(x = mat,
         n = 50,
         fm = "pa",
         rotate = "varimax")

## Summarize heuristic results ----
bind_cols(fa.ns$Analysis,
          fa.nf$vss.stats)

## Set number of factors to extract ----
nFac = 10

## Perform exploratory factor analysis ----
fa.varimax <- fa(r=mat_cor,
                 nfactors = nFac,
                 fm="pa",
                 max.iter=100,
                 rotate="varimax")

# Item Uniqueness: uniqueness, sometimes referred to as noise, corresponds to the proportion of variability, which can not be explained by a linear combination of the factors. A high uniqueness for a variable indicates that the factors do not account well for its variance. (what exists in the unique variance instead of common variance)
fa.varimax.uniqueness = data.table(question = names(fa.varimax$uniquenesses),
                                   uniqueness = fa.varimax$uniquenesses)


# Factor Loadings
fa.varimax.str = fa.varimax$loadings
fa.varimax.pat = promax(fa.varimax$loadings)$loadings

fa.varimax.loadings = list()
for (i in 1:nFac){
  load_pat = fa.varimax.pat[abs(fa.varimax.str[,i]) > 0.01,i]
  load_str = fa.varimax.str[abs(fa.varimax.str[,i]) > 0.01,i]
  fa.varimax.loadings[[i]] = merge(data.table(factor = colnames(fa.varimax.str)[i],
                                              item = names(load_pat),
                                              pattern = load_pat),
                                   data.table(factor = colnames(fa.varimax.str)[i],
                                              item = names(load_str),
                                              structure = load_str),
                                   by = c("factor","item"),
                                   all = T)
}

fa.varimax.loadings = rbindlist(fa.varimax.loadings)

fa.varimax.loadings = fa.varimax.loadings %>% 
  arrange(factor,-abs(structure))

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


for (i in colnames(fa.varimax.scores.parsed)[-1]){
  fa.phe = fa.varimax.scores.parsed %>%
    mutate(phe = get(i)) %>%
    separate(col = dog, 
             into = c("cohort","dog"),
             sep = "\\|") %>%
    mutate(fid = if_else(cohort == "Darwin's Ark",
                         dog,
                         "0"),
           iid = dog) %>%
    select(fid,iid,phe)
  
  write_tsv(fa.phe,
            paste("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/2023-01-01/cross-cohort.DORA.CCDR.",i,".tsv", sep = ""),
            col_names = F)
}


####

geno = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/dat/pca/DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_indep-pairwise_kb-250_r2-0.2.eigenvec")

# Get ages
load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_HLES_dog_owner_v1.0.RData")
HLES_dog_owner$dd_age_years

age_data = HLES_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id, dd_birth_year)

# Canine Eating Behavior (CEB) ~ DORA
load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_AFUS_dog_owner_v1.0.RData")
colnames(AFUS_dog_owner)

answers = AFUS_dog_owner %>%
  pivot_longer(cols = -c("dog_id","owner_id"),
               names_to = "item",
               values_to = "selection",
               values_transform = list(selection = as.character))

colnames(AFUS_dog_owner)[colnames(AFUS_dog_owner) %like% "activities"]

answers %>%
  filter(item %like% "activities") %>%
  filter(dog_id %in% geno$IID) %>%
  ggplot(aes(x = selection)) +
  geom_bar() +
  facet_wrap(.~item, scales = "free")

dap.phe.afus.activity.hunting = AFUS_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id,afus_dd_activities_hunting) %>%
  mutate(phe = if_else(is.na(afus_dd_activities_hunting),
                       0,
                       1)) %>%
  mutate(fid=0,
         iid=dog_id) %>%
  select(fid,iid,phe)
write_tsv(dap.phe.afus.activity.hunting,
          file = "~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/dap.phe.afus.activity.hunting.tsv",
          col_names = F)

dap.phe.afus.activity.field_trial = AFUS_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id,afus_dd_activities_field_trials) %>%
  mutate(phe = if_else(is.na(afus_dd_activities_field_trials),
                       0,
                       1)) %>%
  mutate(fid=0,
         iid=dog_id) %>%
  select(fid,iid,phe)
write_tsv(dap.phe.afus.activity.field_trial,
          file = "~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/dap.phe.afus.activity.field_trial.tsv",
          col_names = F)

dap.phe.afus.activity.agility = AFUS_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id,afus_dd_activities_agility) %>%
  mutate(phe = if_else(is.na(afus_dd_activities_agility),
                       0,
                       1)) %>%
  mutate(fid=0,
         iid=dog_id) %>%
  select(fid,iid,phe)
write_tsv(dap.phe.afus.activity.agility,
          file = "~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/dap.phe.afus.activity.agility.tsv",
          col_names = F)

dap.phe.afus.mdors.costs_money = AFUS_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id,afus_mdors_27_costs_money) %>%
  mutate(phe = afus_mdors_27_costs_money) %>%
  mutate(fid=0,
         iid=dog_id) %>%
  select(fid,iid,phe)
write_tsv(dap.phe.afus.mdors.costs_money,
          file = "~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/dap.phe.afus.mdors.costs_money.tsv",
          col_names = F)

dap.phe.afus.mdors.perceived_cost = AFUS_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id,afus_mdors_subscale_perceived_costs) %>%
  mutate(phe = afus_mdors_subscale_perceived_costs) %>%
  mutate(fid=0,
         iid=dog_id) %>%
  select(fid,iid,phe)
write_tsv(dap.phe.afus.mdors.perceived_cost,
          file = "~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/dap.phe.afus.mdors.perceived_cost.tsv",
          col_names = F)

for (item in colnames(AFUS_dog_owner)[colnames(AFUS_dog_owner) %like% "afus_dora"]){
  write_tsv((AFUS_dog_owner %>%
    filter(dog_id %in% geno$IID) %>%
    mutate(fid=0,
           iid=dog_id,
           phe=get(item)) %>%
    select(fid,iid,phe)),
    paste("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/",
          "dap.phe.",
          item,
          ".tsv",
          sep = ""),
    col_names = F)
}

# Canine Social and Learned Behaviors (CSLB) ~ CCDR
load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_CSLB_v1.0.RData")
geno = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/dat/pca/DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_indep-pairwise_kb-250_r2-0.2.eigenvec")


for (item in colnames(CSLB)[-(1:3)]){
  write_tsv((CSLB %>%
               filter(cslb_year==2020) %>%
               filter(dog_id %in% geno$IID) %>%
               mutate(fid=0,
                      iid=dog_id,
                      phe=get(item)) %>%
               select(fid,iid,phe)),
            paste("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/",
                  "dap.phe.",
                  item,
                  ".tsv",
                  sep = ""),
            col_names = F)
}

library(anytime)
library(lubridate)
age_data = HLES_dog_owner %>%
  filter(dog_id %in% geno$IID) %>%
  select(dog_id, dd_birth_year) %>%
  merge((AFUS_dog_owner %>% 
           select(dog_id,afus_complete_date)), 
        by = "dog_id") %>%
  mutate(birth_year = as.integer(dd_birth_year),
         response_year = year(anydate(afus_complete_date))) %>%
  mutate(age = response_year - birth_year) %>%
  mutate(fid=0,
         iid=dog_id) %>%
  select(fid,iid,age)
write_tsv(age_data,
          "~/Dropbox (UMass Medical School)/Projects/DogAgingProject/gwas/phe/age.tsv", col_names = F)


# Dog Aging Project - Reported Breeds ----
dap.dog = read_csv("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DogAgingMain-BreedsexweightCompar_DATA_LABELS_2022-12-14_1252.csv")
dap.gen = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/dat/pca/DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_indep-pairwise_kb-250_r2-0.2.eigenvec")

dap.seq.breed= dap.dog %>%
  select(dog=`Study ID`,
         breed=`What breed is your dog`) %>%
  filter(dog %in% dap.gen$IID) %>%
  filter(!is.na(breed)) %>%
  group_by(breed) %>%
  summarize(n = n()) %>%
  mutate(pct = n/sum(n))

dap.all.breed= HLES_dog_owner %>%
  mutate(breed=as_character(dd_breed_pure),
         weight=conv_unit(dd_weight_lbs,"lbs","kg")) %>%
  filter(!is.na(breed)) %>%
  select(dog=dog_id,
         breed,
         weight) %>%
  group_by(breed) %>%
  summarize(n = n(),
            weight=mean(weight)) %>%
  ungroup() %>%
  mutate(pct = n/sum(n))


dap.breed = merge(dap.seq.breed,dap.all.breed, by = "breed")
t.test(dap.breed$pct.x,dap.breed$pct.y,paired=T)
cor.test(dap.breed$pct.x,dap.breed$pct.y)

t.test((dap.breed %>% filter(weight > 40) %>% pull(pct.x)),
       (dap.breed %>% filter(weight > 40) %>% pull(pct.y)))

t.test((dap.breed %>% filter(weight < 10) %>% pull(pct.x)),
       (dap.breed %>% filter(weight < 10) %>% pull(pct.y)))

t.test((dap.breed %>% filter(weight < 10) %>% pull(pct.x)),
       (dap.breed %>% filter(weight < 10) %>% pull(pct.y)))

breed.weights = HLES_dog_owner %>%
  mutate(breed=as_character(dd_breed_pure),
         weight=conv_unit(dd_weight_lbs,"lbs","kg")) %>%
  select(dog=dog_id,
         breed,
         weight) %>%
  group_by(breed) %>%
  summarize(breed.weight=mean(weight))

dap.seq.breed= dap.dog %>%
  select(dog=`Study ID`,
         breed=`What breed is your dog`) %>%
  filter(dog %in% dap.gen$IID) %>%
  merge(breed.weights, by = "breed") %>%
  mutate(ext = breed.weight > 40 | breed.weight < 10)

dap.all.breed= HLES_dog_owner %>%
  mutate(breed=as_character(dd_breed_pure),
         weight=conv_unit(dd_weight_lbs,"lbs","kg")) %>%
  select(dog=dog_id,
         breed,
         weight) %>%
  merge(breed.weights, by = "breed") %>%
  mutate(ext = breed.weight > 40 | breed.weight < 10)

t.test(dap.all.breed$ext, dap.seq.breed$ext)

dap.seq.breed= dap.dog %>%
  select(dog=`Study ID`,
         breed=`What breed is your dog`) %>%
  filter(dog %in% dap.gen$IID) %>%
  merge(breed.weights, by = "breed") %>%
  mutate(ext = breed.weight > 40)

dap.all.breed= HLES_dog_owner %>%
  mutate(breed=as_character(dd_breed_pure),
         weight=conv_unit(dd_weight_lbs,"lbs","kg")) %>%
  select(dog=dog_id,
         breed,
         weight) %>%
  merge(breed.weights, by = "breed") %>%
  mutate(ext = breed.weight > 40)

t.test(dap.all.breed$ext, dap.seq.breed$ext)

dap.seq.breed= dap.dog %>%
  select(dog=`Study ID`,
         breed=`What breed is your dog`) %>%
  filter(dog %in% dap.gen$IID) %>%
  merge(breed.weights, by = "breed") %>%
  mutate(ext = breed.weight < 10)

dap.all.breed= HLES_dog_owner %>%
  mutate(breed=as_character(dd_breed_pure),
         weight=conv_unit(dd_weight_lbs,"lbs","kg")) %>%
  select(dog=dog_id,
         breed,
         weight) %>%
  merge(breed.weights, by = "breed") %>%
  mutate(ext = breed.weight < 10)

t.test(dap.all.breed$ext, dap.seq.breed$ext)




dap.seq.breed= dap.dog %>%
  select(dog=`Study ID`,
         breed=`What breed is your dog`) %>%
  filter(dog %in% dap.gen$IID) %>%
  group_by(breed) %>%
  summarize(n = n()) %>%
  mutate(pct = n/sum(n))

dap.all.breed= HLES_dog_owner %>%
  mutate(breed=as_character(dd_breed_pure),
         weight=conv_unit(dd_weight_lbs,"lbs","kg")) %>%
  select(dog=dog_id,
         breed,
         weight) %>%
  group_by(breed) %>%
  summarize(n = n(),
            weight=mean(weight)) %>%
  ungroup() %>%
  mutate(pct = n/sum(n))

p = merge(dap.seq.breed,dap.all.breed, by = "breed") %>%
  merge(breed.weights, by = "breed") %>%
  mutate(strata = cut(breed.weight,
                      breaks= c(0,10,20,30,40,max(breed.weight)),
                      include.lowest = T)) %>%
  ggplot(aes(x = pct.x,
             y = pct.y,
             label = breed)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(color = strata,
                 shape = strata)) +
  geom_text_repel()+
  labs(x = "proportion of sequenced dogs",
       y = "proportion of surveyed dogs") +
  scale_x_log10()+
  scale_y_log10()+
  theme_pubr()

ggsave(filename="~/Dropbox (UMass Medical School)/Thesis/figures/DAP_BREED_PROP.pdf",
       plot = p,
       width =6*2,
       height= 6*2,
       device = "pdf")
library(ggrepel)


da.breeds = dogs %>%
  filter(dog %in% darwinsark.geno$dog) %>%
  group_by(breed) %>%
  summarize(n= n()) %>%
  ungroup() %>%
  mutate(pct = n/sum(n))

dap.breeds = dap.dog %>% 
  filter(`Study ID` %in% dap.gen$IID) %>%
  mutate(breed = tolower(`What breed is your dog`)) %>%
  select(breed) %>%
  group_by(breed) %>%
  summarize(n= n()) %>%
  ungroup() %>%
  mutate(pct = n/sum(n))

breeds = merge(da.breeds,dap.breeds,by="breed") %>% filter(!is.na(breed))
cor.test(breeds$pct.x,breeds$pct.y)
cor.test(breeds$pct.x,breeds$pct.y)$p.value

library(corrplot)

breed.mat = breeds %>% select(pct.x,pct.y)
rownames(breed.mat) = (breeds %>% filter())
cor((breeds %>% select(pct.x,pct.y)))

p=breeds %>%
  merge((breed.weights %>% mutate(breed= tolower(breed))), by = "breed") %>%
  mutate(strata = cut(breed.weight,
                      breaks= c(0,
                                10,
                                20,
                                30,
                                40,
                                max(breed.weight)),
                      include.lowest = T)) %>%
  filter(!is.na(breed)) %>%
  ggplot(aes(x = pct.x,
             y = pct.y,
             label = breed)) +
  geom_point(aes(color = strata,
                 shape = strata)) +
  geom_text_repel() +
  scale_shape_manual(values = c(6, 5, 0, 12, 2)) +
  scale_x_log10()+
  scale_y_log10()+
  labs(x = "Darwin's Ark\nproportion of breeds reported",
       y = "Dog Aging Project\nproportion of breeds reported") +
  theme_pubr()

ggsave(filename="~/Dropbox (UMass Medical School)/Thesis/figures/DAP_VS_DA_BREED_PROP.pdf",
       plot = p,
       width =6,
       height= 6,
       device = "pdf")

ggplotly(p)


p=da.breeds %>%
  mutate(cohort= "Darwin's Ark",
         pct = -pct) %>%
  filter(breed %in% dap.breeds$breed) %>%
  bind_rows((dap.breeds %>%
               filter(breed %in% da.breeds$breed) %>%
               mutate(cohort="Dog Aging Project"))) %>%
  merge((breed.weights %>% 
           mutate(breed= tolower(breed))), by = "breed") %>%
  mutate(strata = cut(breed.weight,
                      breaks= c(0,
                                10,
                                20,
                                30,
                                40,
                                max(breed.weight)),
                      include.lowest = T)) %>%
  filter(!is.na(breed)) %>%
  ggplot(aes(x = reorder(breed,pct),
             y = pct,
             fill = cohort)) +
  coord_flip() +
  geom_col() +
  theme_pubr()

da.breeds %>%
  mutate(q25 = quantile(pct,.5)) %>%
  filter(pct<q25)

breeds %>%
  mutate(fold = pct.y/pct.x) %>%
  merge((breed.weights %>% mutate(breed= tolower(breed))), by = "breed") %>% View()
  mutate(strata = cut(breed.weight,
                      breaks= c(0,
                                10,
                                20,
                                30,
                                40,
                                max(breed.weight)),
                      include.lowest = T)) %>%
  filter(n.x > 10 | n.y > 10) %>%
  filter(!is.na(breed)) %>%
  ggplot(aes(x = reorder(breed,fold),
             y = fold,
             fill = strata)) +
  coord_flip() +
  geom_col() +
  labs(x = "breed",
       y = "proportion in Dog Aging Project /\nproportion in Darwin's Ark") +
  theme_pubr()
