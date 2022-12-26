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

args <- parser$parse_args()
print(args)

### Set working directory ----
workDir = args$dir
setwd(workDir)

### Set data freeze ----
freezeDate = args$date
freezeDir = paste(workDir,"data","database","data_freeze",freezeDate, sep = "/")

### Load from saved image, if available ----
load(paste(workDir, "data", paste(paste("DarwinsArk", freezeDate, sep="_"), ".RData", sep=""), sep="/"))

# LOAD AND PARSE DATA ----

## DOGS ----
dogs_raw = read_csv(paste(freezeDir,
                          "dogs.csv",
                          sep = "/"),
                    col_names = TRUE,
                    na = c("NA",
                           "NULL",
                           "null",
                           "")) %>% as.data.table()

### Protocol for data wrangling ----
# 1. Remove fake dogs enrolled for testing
# 2. Assign columns more meaningful names
# 3. Select columns in desired order

### Select and parse valid entries from `dogs` table ----
dogs = dogs_raw[fake == 0 & !name %like% "Test",
                list(dog = id,
                     owner,
                     group = group_id,
                     sex,
                     neutered,
                     name,
                     text_age = age,
                     birth_date = anydate(birthday),
                     consent_date = anydate(consent_date),
                     flagged_deceased_date = anydate(as.numeric(deceased)), 
                     breed1,
                     breed1_id,
                     breed2,
                     breed2_id,
                     purebred,
                     ancestry_owner_feedback_agree = is_breed_mix_right,
                     ancestry_owner_feedback_date = feedback_submit_date,
                     ancestry_owner_feedback_known = is_breeds_known,
                     ancestry_owner_feedback_source = breed_info_source,
                     ancestry_owner_feedback_source_note = breed_info_source_others,
                     ancestry_owner_feedback_info = ancestry_info,
                     profile_public = is_public,
                     image_id = image,
                     image_private = image_private,
                     working_dog = is_working
                )]

## FORMATS ----
formats = read_csv(paste(freezeDir,
                         "formats.csv",
                         sep = "/"),
                   col_names = TRUE,
                   na = c("NA",
                          "NULL",
                          "null",
                          "")) %>% as.data.table()

formats.options = list()
i = 1
for (j in formats$id){
  dt = data.table(id = j,
                  index = as.character(seq(0,length(unlist(strsplit(formats[id == j]$options, "|", fixed = T)))-1)),
                  option = unlist(strsplit(formats[id == j]$options, "|", fixed = T)))
  formats.options[[i]] = dt
  i = i + 1
}
formats.options = rbindlist(formats.options) %>% as.data.table() %>% na.omit()

## SURVEYS ----
surveys = read_csv(paste(freezeDir,
                         "surveys.csv",
                         sep = "/"),
                   col_names = TRUE,
                   na = c("NA",
                          "NULL",
                          "null",
                          "")) %>% as.data.table()

## QUESTIONS ----
questions = read_csv(paste(freezeDir,
                           "questions.csv",
                           sep = "/"),
                     col_names = TRUE,
                     na = c("NA",
                            "NULL",
                            "null",
                            "")) %>% 
  merge(formats,
        by.x = "format",
        by.y = "id",
        all.x = T) %>%
  merge(surveys,
        by.x = "survey",
        by.y = "id",
        all.x = T) %>%
  as.data.table()

## ANSWERS ----
answers_raw = read_csv(file = paste(freezeDir,
                                    "answers.csv",
                                    sep = "/"),
                       col_names = TRUE,
                       na = c("NA",
                              "NULL",
                              "null")) %>% as.data.table()

### Protocol for data wrangling ----
# 1. Remove fake dogs enrolled for testing
# 2. Convert empty entries into missing (NA)
# 3. Assign `format` keyed on `question`
# 4. For multiple responses, select latest time stamp
# 5. For multiple choice responses, split into rows
# 6. Assign `option` from `answer` values, where possible

### Select and parse valid entries from `answers` ----
answers = answers_raw %>%
  filter(dog %in% dogs$dog) %>%
  mutate(answer = if_else(answer == "",
                          NA_character_,
                          answer)) %>%
  filter(!is.na(answer)) %>%
  merge((questions %>% select(id,format)),
        by.x = "question",
        by.y = "id",
        all.x = T) %>%
  mutate(multi = question %in% (answers_raw %>% 
                                  filter(grepl("\\|",answer)) %>% 
                                  arrange(question) %>%
                                  pull(question) %>%
                                  unique())) %>%
  group_by(dog,user,question) %>%
  filter(time == max(time)) %>%
  slice_tail(n = 1) %>% 
  separate_rows(answer,
                sep = "\\|",
                convert = T) %>%
  group_by(dog,user,question) %>%
  mutate(nest = row_number()-1) %>%
  as.data.table()

answers[!question %in% (questions %>% 
                          filter(format %in% c(3,4)) %>% 
                          pull(id)), 
        index := as.character(answer)]

answers[question %in% c(111,112), 
        index := as.character(nest)]

### Assign options ----
answers = answers %>%
  merge(formats.options,
        by.x = c("format","index"),
        by.y = c("id","index"),
        all.x = T)

### Remove duplicate responses ----
answers = answers %>%
  group_by(dog,question) %>%
  mutate(n = n()) %>% 
  mutate(row = row_number()) %>%
  mutate(nopt = max(nest)==0) %>%
  filter((nopt & row == max(row) & n>1) | (!nopt) | (nopt & n == 1)) %>%
  select(-n,-nest,-nopt,-row)

### Sort answers table ----
answers = answers %>%
  arrange(dog,
          user,
          question) %>%
  as.data.table()

## REGIONS ----

### Read raw table `zipcodes` ----
zips = read_csv(file = paste(freezeDir,"zipcodes.csv", sep = "/"),
                na = c("NA",
                       "NULL", 
                       "null",
                       "")) %>% as.data.table()

### Unify state and country ----
locs = zips %>%
  select(owner = user_id,
         raw_state=state, shipping_state, billing_state,
         raw_country=country, shipping_country, billing_country) %>%
  pivot_longer(cols = c("raw_state",
                        "shipping_state",
                        "billing_state",
                        "raw_country",
                        "shipping_country",
                        "billing_country"),
               names_to = "source",
               values_to = "value") %>%
  mutate(type = if_else(source %like% "state",
                        "state",
                        "country")) %>%
  group_by(owner,type) %>%
  filter(!is.na(value)) %>%
  mutate(consensus = length(unique(value)) == 1)

locs = bind_rows((locs %>%
                    filter(consensus == TRUE) %>%
                    select(owner,type,value) %>%
                    unique() %>%
                    pivot_wider(id_cols = "owner", 
                                names_from = "type",
                                values_from = "value") %>%
                    mutate(country = if_else(country %in% c("United States","US"),
                                             "USA",
                                             country))),
                 (locs %>%
                    filter(source %in% c("shipping_state","billing_country") & consensus == FALSE) %>%
                    pivot_wider(id_cols = "owner", 
                                names_from = "type",
                                values_from = "value") %>%
                    mutate(country = if_else(country %in% c("United States","US"),
                                             "USA",
                                             country)))) %>%
  pivot_longer(cols = c("state","country"),
               names_to = "source",
               values_to = "value") %>%
  filter(!is.na(value)) %>% 
  pivot_wider(id_cols = "owner",
              names_from = "source",
              values_from = "value")

### If state matches and country empty, fill USA ----
locs = locs %>%
  mutate(country = if_else(is.na(country) & !is.na(state) & state %in% state.abb,
                           "USA",
                           country))

### Unify zips ----
zips = zips %>% select(owner = user_id,
                       raw_zip = postcode,
                       bill_zip = billing_postcode,
                       ship_zip = shipping_postcode) %>%
  pivot_longer(cols = c("raw_zip",
                        "bill_zip",
                        "ship_zip"),
               names_to = "source",
               values_to = "value") %>%
  group_by(owner) %>%
  filter(!is.na(value)) %>%
  mutate(consensus = length(unique(value)) == 1)

zips = bind_rows((zips %>%
                    filter(consensus == TRUE) %>%
                    select(owner,zip=value) %>%
                    unique()),
                 (zips %>%
                    filter(consensus == FALSE & source == "ship_zip") %>%
                    select(owner,zip=value))) %>%
  as.data.table()

### Repair zips in locs ----
locs = locs %>% merge(zips, by = "owner", all = T)

### Remove zips from local env ----
rm(zips)

### Match to dogs ----
locs = dogs %>%
  select(owner,dog) %>%
  merge(locs, by = "owner") %>%
  select(-owner)

### If responses to Q#252 "Enter the country and postal code of your home (where DOG spends the most time):", then obtain and compare zip codes ----
locs = locs %>%
  bind_rows(answers %>% 
              filter(question == 252) %>%
              select(dog,answer) %>%
              mutate(parsed = map(answer, fromJSON)) %>%
              unnest_wider(parsed) %>%
              select(dog,zip,country,other) %>%
              mutate(country = if_else(country == "United States of America",
                                       "USA",
                                       country))) %>%
  select(-other) %>%
  as.data.table()

### Get state from zip ----
locs = locs %>%
  mutate(zip = if_else(country == "USA",
                       normalize_zip(zip),
                       zip)) %>%
  mutate(state = if_else(!is.na(zip) & is.na(state) & country == "USA",
                         reverse_zipcode(zip)$state,
                         state)) %>%
  as.data.table()

### Collapse dogs ----
locs = locs %>%
  group_by(dog) %>%
  summarize(state = unique(state[!is.na(state)])[1],
            country = unique(country[!is.na(country)])[1],
            zip = unique(zip[!is.na(zip)])[1]) %>% 
  as.data.table()

### Define regions of United States ----
locs[state == 'AK', region := "Noncontiguous (Alaska)"]
locs[state == 'HI', region := "Noncontiguous (Hawaii)"]
locs[state == 'PR', region := "Noncontiguous (Puerto Rico)"]
locs[state == 'AS', region := "Noncontiguous (American Samoa)"]
locs[state %in% c("CT", "ME", "MA", "NH", "RI", "VT"), region := "Northeast (New England)"]
locs[state %in% c("NJ", "NY", "PA"), region := "Northeast (Mid-Atlantic)"]
locs[state %in% c("DE", "FL", "GA", "MD", "NC", "SC", "VT", "DC", "WV", "VA"), region := "South (South Atlantic)"]
locs[state %in% c("AL", "KY", "MS", "TN"), region := "South (East South Central)"]
locs[state %in% c("AR", "LA", "OK", "TX"), region := "South (West South Central)"]
locs[state %in% c("IL", "IN", "MI", "OH", "WI"), region := "Midwest (East North Central)"]
locs[state %in% c("IA", "KS", "MN", "MO", "NE", "ND", "SD"), region := "Midwest (West North Central)"]
locs[state %in% c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY"), region := "West (Mountain)"]
locs[state %in% c("CA", "OR", "WA"), region := "West (Pacific)"]

### Get census data ----
census = read_csv(file = paste(workDir,
                               "/ref/",
                               "US_census_2010.csv", 
                               sep="")) %>% 
  as.data.table()
colnames(census) = c("zip",
                     "population",
                     "landmass",
                     "density")

### Define environ based on densities ----
locs[zip %in% (census %>% filter(density >= 700) %>% pull(zip)), environ := "urban"]
locs[zip %in% (census %>% filter(density >= 100 & density < 700) %>% pull(zip)), environ := "suburban"]
locs[zip %in% (census %>% filter(density < 100) %>% pull(zip)), environ := "rural"]

### Define environ based on Q#254 ----
locs = locs %>% 
  merge((answers %>%
           filter(question == 254) %>%
           select(dog,answer)),
        by = "dog",
        all.x = T) %>% as.data.table()

locs[answer==0, environ := "urban"]
locs[answer==1, environ := "suburban"]
locs[answer==2, environ := "rural"]

### Assign region and environ data ----
dogs = dogs %>%
  merge((locs %>% 
           select(dog,
                  state,
                  country,
                  region,
                  environ)),
        by = "dog",
        all.x = T) %>% 
  arrange(dog) %>%
  as.data.table()

## ORIGINS ----

### Assign `adopt_date` to `dogs` ---
dogs = dogs %>%
  merge((answers %>%
           filter(question == 116) %>%
           select(dog,answer) %>%
           mutate(adopt_date = anydate(answer)) %>%
           select(dog,adopt_date)),
        by = "dog",
        all.x = T)

### Assign `adopt_age` to `dogs` ----
dogs = dogs %>%
  merge((answers %>%
           filter(question == 241) %>%
           select(dog,adopt_age=option)),
        by = "dog",
        all.x = T)

### Assign `origin` to `dogs` ----
dogs = dogs %>%
  merge((answers %>%
           filter(question == 117) %>%
           select(dog,origin=option)))

## DATES AND TIMES ----
### Select start as earliest response time stamp ----
dogs = dogs %>%
  merge((answers_raw %>%
           group_by(dog) %>%
           arrange(-desc(time)) %>%
           slice(1) %>%
           select(dog,time) %>%
           mutate(start_date = anydate(time)) %>%
           select(dog,start_date)),
        by = "dog",
        all.x = T)

### Assign year events ----
dogs = dogs %>%
  mutate(birth_year = year(birth_date),
         adopt_year = year(adopt_date),
         consent_year = year(consent_date),
         start_year = year(start_date),
         death_year = year(flagged_deceased_date)
  ) %>% 
  arrange(dog) %>%
  as.data.table()

### If adopted as young puppy, set born same year ----
dogs[is.na(birth_year) & !is.na(adopt_year) & adopt_age %like% "Younger puppy",
     birth_year := adopt_year]

### If adopted as older puppy, set born -1 year ----
dogs[is.na(birth_year) & !is.na(adopt_year) & adopt_age %like% "Older puppy",
     birth_year := adopt_year-1]

### If adopted as young adult, set born -2 years ----
dogs[is.na(birth_year) & !is.na(adopt_year) & adopt_age %like% "Young adult",
     birth_year := adopt_year-2]

### Remove illogical dates and years ----
dogs[birth_year < 1980, birth_date := NA]
dogs[birth_year < 1980, birth_year := NA]
dogs[birth_year > year(anydate(freezeDate)), birth_date := NA]
dogs[birth_year > year(anydate(freezeDate)), birth_year := NA]

dogs[adopt_year < 1980, adopt_date := NA]
dogs[adopt_year < 1980, adopt_year := NA]
dogs[adopt_year > year(anydate(freezeDate)), adopt_date := NA]
dogs[adopt_year > year(anydate(freezeDate)), adopt_year := NA]

dogs[start_year < 1980, start_date := NA]
dogs[start_year < 1980, start_year := NA]
dogs[start_year > year(anydate(freezeDate)), start_date := NA]
dogs[start_year > year(anydate(freezeDate)), start_year := NA]

## AGES ----

### In `answers`, calculate age at response ----
answers = answers %>%
  merge((dogs %>%
           select(dog,
                  birth_year)),
        by = "dog") %>%
  mutate(age = year(anydate(time)) - birth_year) %>%
  mutate(age = replace(age, age<=0, NA)) %>%
  as.data.table()

### In `dogs`, assign average age and age group ----
dogs = dogs %>% 
  merge(
    (answers %>%
       group_by(dog) %>%
       summarize(age = ceiling(mean(age, na.rm = T)),
                 age_group = if_else(mean(age, na.rm = T) >= 1,
                                     "adult",
                                     "puppy",
                                     "adult")) %>%
       filter(!is.nan(age))),
    all.x = T,
    by = "dog") %>%
  as.data.table()

## SEX AND SPAY/NEUTER ----
dogs[, sex := tolower(sex)]
dogs[, neutered := tolower(neutered)]
dogs[neutered != "yes" & neutered != "no", neutered := NA]

dogs[sex == "male" & neutered == "yes", sex_status :=  "neutered male"]
dogs[sex == "male" & neutered == "no", sex_status :=  "intact male"]
dogs[sex == "male" & is.na(neutered), sex_status :=  "unspecified male"]
dogs[sex == "female" & neutered == "yes", sex_status :=  "spayed female"]
dogs[sex == "female" & neutered == "no", sex_status :=  "intact female"]
dogs[sex == "female" & is.na(neutered), sex_status :=  "unspecified female"]

## PEDIGREES ----
# 1 = registered purebred, 0 = not registered, -1 = unknown
dogs[, purebred := tolower(purebred)]
dogs[purebred == "1", purebred := "yes"]
dogs[purebred != "yes" & purebred != "no", purebred := NA]

## BREEDS ----

### Read raw data table `breeds` ----
breeds = read_csv(file = paste(freezeDir,"breeds.csv", sep = "/"),
                  na = c("NA", 
                         "NULL", 
                         "null")) %>% as.data.table()
breeds[, breed_name := tolower(breed_print)]

### 1. Save reported breed entries ----
dogs[, breed1.original := breed1]
dogs[, breed2.original := breed2]

### 2. Synchronize strings ----
dogs[, breed1 := tolower(iconv(enc2utf8(breed1), sub = "byte"))]
dogs[, breed2 := tolower(iconv(enc2utf8(breed2), sub = "byte"))]

dogs[breed1 %in% c("","?","null","NULL"), 
     breed1 := NA]
dogs[breed2 %in% c("","?","null","NULL"), 
     breed2 := NA]

### 3. Handle `breed1_id` and `breed2_id` ----
dogs[breed1_id == 0, breed1_id := NA]
dogs[breed2_id == 0, breed2_id := NA]

### 4. Match IDs to breeds ----
breed1_id_match = merge(dogs[is.na(breed1) & !is.na(breed1_id)], 
                        breeds, 
                        by.x = "breed1_id", 
                        by.y = "id", 
                        all.x = T)[, 
                                   c("dog",
                                     "breed1_id",
                                     "breed_name"
                                   )]
setkey(breed1_id_match,dog)
setkey(dogs,dog)
dogs[is.na(breed1) & !is.na(breed1_id)]$breed1 <- breed1_id_match$breed_name

breed2_id_match = merge(dogs[is.na(breed2) & !is.na(breed2_id)], 
                        breeds, 
                        by.x = "breed2_id", 
                        by.y = "id", 
                        all.x = T)[, 
                                   c("dog",
                                     "breed2_id",
                                     "breed_name"
                                   )]
setkey(breed2_id_match,dog)
setkey(dogs,dog)
dogs[is.na(breed2) & !is.na(breed2_id)]$breed2 <- breed2_id_match$breed_name

### 5. Drop `breed1_id` and `breed2_id` ----
dogs[, breed1_id := NULL]
dogs[, breed2_id := NULL]

### 6. Move `breed2` to `breed1` if `breed1` is NA ----
replace_dogs = dogs[is.na(breed1) & !is.na(breed2)]$id
dogs[dog %in% replace_dogs, breed1 := breed2]
dogs[dog %in% replace_dogs, breed2 := NA]
rm(replace_dogs)

### 7. If `breed1` == `breed2`, then drop `breed2` ----
dogs[breed1==breed2, breed2 := NA]

### 8. Set alternate breed conversions ----
breed_alternates = c(
  "goldendoodle" = "goldendoodle (golden retriever x poodle)",
  "golden doodle" = "goldendoodle (golden retriever x poodle)",
  "golden-doodle" = "goldendoodle (golden retriever x poodle)",
  "labradoodle" = "labradoodle (labrador retriever x poodle)",
  "labra doodle" = "labradoodle (labrador retriever x poodle)",
  "labra-doodle" = "labradoodle (labrador retriever x poodle)",
  "aussiedoodle" = "aussiedoodle (australian shepherd x poodle)",
  "schnoodle" = "schnoodle (schnauzer x poodle)",
  "cockapoo" = "cockapoo (cocker spaniel x poodle)",
  "cock-a-poo" = "cockapoo (cocker spaniel x poodle)",
  "cocka-poo" = "cockapoo (cocker spaniel x poodle)",
  "yorkiepoo" = "yorkiepoo (yorkshire terrier x poodle)",
  "yorkie-poo" = "yorkiepoo (yorkshire terrier x poodle)",
  "maltipoo" = "maltipoo (maltese x poodle)",
  "malti-poo" = "maltipoo (maltese x poodle)",
  "peek-a-poo" = "pekapoo (pekingese x poodle)",
  "pekeapoo" = "pekapoo (pekingese x poodle)",
  "pekapoo" = "pekapoo (pekingese x poodle)",
  "shih-poo" = "shih-poo (shih tzu x poodle)",
  "bichpoo" = "poochon (bichon x poodle)",
  "longhaired dachshund" = "longhaired dachshund",
  "wirehaired dachshund" = "wirehaired dachshund",
  "smooth dachshund" = "smooth dachshund",
  "chihuahua (long coat)" = "longhaired chihuahua",
  "chihuahua (long" = "longhaired chihuahua",
  "long hair chihuahua" = "longhaired chihuahua",
  "miniature australian shepherd" = "miniature american shepherd",
  "miniature american shepherd" = "miniature american shepherd",
  "american shepherd" = "miniature american shepherd",
  "miniature bull terrier" = "miniature bull terrier",
  "miniature poodle" = "miniature poodle",
  "moyen poodle" = "miniature poodle",
  "toy poodle" = "toy poodle",
  "standard poodle" = "standard poodle",
  "rough collie" = "rough collie",
  "smooth collie" = "smooth collie",
  "parson russell terrier" = "parson russell terrier",
  "russell terrier" = "russell terrier",
  "llewellyn setter" = "llewellin setter",
  "llewellyn" = "llewellin setter",
  "llewellin setter" = "llewellin setter",
  "llewellin" = "llewellin setter",
  "belgian malinois" = "belgian malinois",
  "belgian groenendael" = "belgian groenendael",
  "belgian sheepdog" = "belgian groenendael",
  "belgian laekenois" = "belgian laekenois",
  "belgian tervuren" = "belgian tervuren",
  "belgian shepherd" = "belgian tervuren",
  "belgian shepherd tervuren" = "belgian tervuren",
  "malinois" = "belgian malinois",
  "tervuren" = "belgian tervuren",
  "laekenois" = "belgian laekenois",
  "groenendael" = "belgian groenendael",
  "phalene" = "phalene",
  "papillon" = "papillon"
  
)

### 9. Assign `breed.alternate` ----
dogs$breed.alternate <- ifelse(
  str_detect(dogs$breed1,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_alternates)), collapse = "|")), 
  breed_alternates[str_match(dogs$breed1, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_alternates)), collapse = "|"))[,1]],
  ifelse(
    str_detect(dogs$breed2,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_alternates)), collapse = "|")), 
    breed_alternates[str_match(dogs$breed2, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_alternates)), collapse = "|"))[,1]],
    NA_character_))

### 10. Breed string conversions ----
breed_corrections <- c(
  # NA (mixed-breed dog)
  "mix" = NA,
  "mutt" = NA,
  "cross" = NA,
  "%" = NA,
  "50" = NA,
  "25" = NA,
  "+" = NA,
  "?" = NA,
  "&" = NA,
  "half" = NA,
  "maybe" = NA,
  "other" = NA,
  # poodle
  "miniature poodle" = "poodle",
  "moyen poodle" = "poodle",
  "toy poodle" = "poodle",
  "standard poodle" = "poodle",
  "poodle" = "poodle",
  # NA (cross-breed dog)
  "doodle" = NA,
  "oodle" = NA,
  "poo" = NA,
  # dachshund
  "longhaired dachshund" = "dachshund",
  "wirehaired dachshund" = "dachshund",
  "smooth dachshund" = "dachshund",
  "dachshund" = "dachshund",
  # chihuahua
  "chihuahua (long coat)" = "chihuahua",
  "chihuahua (long" = "chihuahua",
  "long hair chihuahua" = "chihuahua",
  "chihuahua" = "chihuahua",
  # australian shepherd
  "aussie" = "australian shepherd",
  "miniature australian shepherd" = "australian shepherd",
  "miniature american shepherd" = "australian shepherd",
  "american shepherd" = "australian shepherd",
  "australian shepherd" = "australian shepherd",
  # collie
  "rough collie" = "collie",
  "smooth collie" = "collie",
  # russell terrier
  "parson russell terrier" = "russell terrier",
  "jack russell terrier" = "russell terrier",
  "russell terrier" = "russell terrier",
  # american pit bull terrier
  "APBT" = "american pit bull terrier",
  "pit bull" = "american pit bull terrier",
  "pitty" = "american pit bull terrier",
  "pitbull terrier" = "american pit bull terrier",
  # american bully
  "american bully" = "american bully",
  # german shepherd dog
  "GSD" = "german shepherd dog",
  "german shepherd" = "german shepherd dog",
  "german shepard" = "german shepherd dog",
  "german sherpherd" = "german shepherd dog",
  "german shepherd dog" = "german shepherd dog",
  "white shepherd" = "german shepherd dog",
  "american white shepherd" = "german shepherd dog",
  "berger blanc suisse shepherd" = "german shepherd dog",
  "berger blanc suisse" = "german shepherd dog",
  "black shepherd" = "german shepherd dog",
  # papillon
  "phalene" = "papillon",
  "papillon" = "papillon",
  # caucasian shepherd dog
  "russian shepherd" = "caucasian shepherd dog",
  "russian shepard" = "caucasian shepherd dog",
  "ovcharka" = "caucasian shepherd dog",
  "caucasian shepherd dog" = "caucasian shepherd dog",
  # dutch shepherd
  "dutch shepherd sheepdog" = "dutch shepherd",
  # belgian shepherd
  "belgian malinois" = "belgian shepherd",
  "belgian groenendael" = "belgian shepherd",
  "belgian sheepdog" = "belgian shepherd",
  "belgian laekenois" = "belgian shepherd",
  "belgian tervuren" = "belgian shepherd",
  "belgian shepherd" = "belgian shepherd",
  "belgian shepherd tervuren" = "belgian shepherd",
  "malinois" = "belgian shepherd",
  "groenendael" = "belgian shepherd",
  "laekenois" = "belgian shepherd",
  "tervuren" = "belgian shepherd",
  # labrador retriever
  "lab" = "labrador retriever",
  "fox red lab" = "labrador retriever",
  "labra" = "labrador retriever",
  "labrador retriever" = "labrador retriever",
  # golden retriever
  "golden" = "golden retriever",
  "golden retriever" = "golden retriever",
  # cocker spaniel
  "english cocker" = "cocker spaniel",
  "american cocker" = "cocker spaniel",
  "cocker" = "cocker spaniel",
  # american eskimo dog
  "american eskimo" = "american eskimo dog",
  # basset hound
  "bassett" = "basset hound",
  "bassett hound" = "basset hound",
  # biewer terrier
  "biewer" = "biewer terrier",
  # bluetick coonhound
  "blue tick hound" = "bluetick coonhound",
  # border collie
  "boarder collie" = "border collie",
  # puli
  "hungarian puli" = "puli",
  # catahoula leopard dog
  "louisiana catahoula leopard dog" = "catahoula leopard dog",
  # black and tan coonhound
  "black & tan coon hound" = "black and tan coonhound",
  "black & tan coonhound" = "black and tan coonhound",
  # american english coonhound
  "red tick coon hound" = "american english coonhound",
  "red tick coonhound" = "american english coonhound
",
  "redtick coonhound" = "american english coonhound",
  # redbone coonhound
  "red bone coon hound" = "redbone coonhound",
  "red bone coonhound" = "redbone coonhound",
  # black mouth cur
  "black mouth curr" = "black mouth cur",
  "black mouth mountain cur" = "black mouth cur",
  "mountain cur" = "black mouth cur",
  "mountain curr" = "black mouth cur",
  # american bulldog
  "amerian bulldog" = "american bulldog",
  # english bulldog
  "olde english bulldogge" = "english bulldog",
  # english setter
  "llewellin" = "english setter",
  "llewellin setter" = "english setter",
  "llewellyn setter" = "english setter",
  "llewellyn setter" = "english setter",
  # brittany
  "brittany spaniel" = "brittany",
  "brittany" = "brittany",
  # bull terrier,
  "miniature bull terrier" = "bull terrier",
  # cane corso
  "cane corso italiano" = "cane corso",
  # bouvier des flandres
  "bouvier" = "bouvier des flandres",
  # australian cattle dog
  "heeler" = "australian cattle dog",
  "cattle dog" = "australian cattle dog",
  # kai ken
  "kai dog" = "kai ken",
  "kai ken" = "kai ken",
  # bolonka
  "russian tsvetnaya bolonka" = "bolonka",
  "bolonka" = "bolonka",
  # lagotto romagnolo
  "lagotti romagnoli" = "lagotto romagnolo",
  "lagotto romagnoli" = "lagotto romagnolo",
  "lagotti romangnoli" = "lagotto romagnolo",
  "lagotto romangnoli" = "lagotto romagnolo",
  # akbash dog
  "akbash" = "akbash dog",
  "akbash dog" = "akbash dog",
  # spinoni italiani
  "spinone" = "spinoni italiani",
  "spinone italiano" = "spinoni italiani",
  # taiwan dog
  "taiwanese (formosan) mountain dog" = "taiwan dog",
  "thaiwaneese mountian dog" = "taiwan dog",
  "formosan mountain dog" = "taiwan dog",
  # greyhound
  "nga greyhound" = "greyhound",
  # wire fox terrier
  "wire hair fox terrier" = "wire fox terrier",
  "wire haired fox terrier" = "wire fox terrier",
  "wirehaired fox terrier" = "wire fox terrier",
  # newfoundland
  "nueflend" = "newfoundland",
  "newfie" = "newfoundland",
  "newfoundland" = "newfoundland",
  # others
  "canaan" = "canaan dog",
  "german wire hair" = "german wirehaired pointer",
  "king charles" = "cavalier king charles spaniel",
  "scottish deerhound" = "scottish deerhound",
  "scottish terrier" = "scottish terrier",
  "shih tzu" = "shih tzu",
  "silky terrier" = "silky terrier",
  "west highland terrier" = "west highland terrier",
  "westie" = "west highland terrier",
  "wheaten terrier" = "soft coated wheaten terrier",
  "yorkie" = "yorkshire terrier",
  "flat-coated retrievers" = "flat-coated retriever",
  "nz hunterway" = "huntaway",
  "patterdale" = "patterdale terrier"
)

### 11. Correct `breeds` table ----
breeds$breed_name_general <- ifelse(
  str_detect(breeds$breed_name,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_corrections)), collapse = "|")), 
  breed_corrections[str_match(breeds$breed_name, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_corrections)), collapse = "|"))[,1]],
  breeds$breed_name)

### 12. Correct `breed1` and `breed2` entries ----
dogs$breed1 <- ifelse(
  str_detect(dogs$breed1,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_corrections)), collapse = "|")), 
  breed_corrections[str_match(dogs$breed1, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_corrections)), collapse = "|"))[,1]],
  dogs$breed1)

dogs$breed2 <- ifelse(
  str_detect(dogs$breed2,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_corrections)), collapse = "|")), 
  breed_corrections[str_match(dogs$breed2, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(breed_corrections)), collapse = "|"))[,1]],
  dogs$breed2)

### 13. Non-breed strings ----
nonbreeds = c(
  "africanis" = "village dog (South Africa)",
  "bush" = "village dog (South Africa)",
  "sato" =      "village dog (Puerto Rico)",
  "soi" =       "village dog (Thailand)",
  "potcake" =   "village dog (Caribbean)",
  "pot cake" =  "village dog (Caribbean)",
  "kampung" =   "village dog (Indonesia)",
  "desi dog" =  "village dog (India)",
  "desi kutta" =  "village dog (India)",
  "indog" = "village dog (India)",
  "Indian pariah" =  "village dog (India)",
  "pariah" =   "non-breed or village dog (any)",
  "village" =   "non-breed or village dog (any)",
  "street" =    "non-breed or village dog (any)",
  "feral" =     "non-breed or village dog (any)",
  "stray" =     "non-breed or village dog (any)",
  "local" =     "non-breed or village dog (any)",
  "native" =    "non-breed or village dog (any)",
  "community" = "non-breed or village dog (any)",
  "free" =      "non-breed or village dog (any)",
  "homeless" =  "non-breed or village dog (any)",
  "island" =  "non-breed or village dog (any)",
  "dingo" =  "non-breed or village dog (any)"
)

### 14. Assign `breed.nonbreed` ----
dogs$breed.nonbreed <- ifelse(
  str_detect(dogs$breed1,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(nonbreeds)), collapse = "|")), 
  nonbreeds[str_match(dogs$breed1, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(nonbreeds)), collapse = "|"))[,1]],
  ifelse(str_detect(dogs$breed2,paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(nonbreeds)), collapse = "|")), 
         nonbreeds[str_match(dogs$breed2, paste(gsub("([^a-zA-Z0-9 ])", "\\\\\\1",names(nonbreeds)), collapse = "|"))[,1]],
         NA_character_))

### 15. Remove `breed1` and `breed2` for nonbreed ----
dogs[!is.na(breed.nonbreed), breed1 := NA]
dogs[!is.na(breed.nonbreed), breed2 := NA]

### 16. Assign `breed` for dogs with only `breed1` ----
dogs[!is.na(breed1) & is.na(breed2), breed := breed1]

### 17. Drop `breed1` and `breed2` ----
dogs[, breed1 := NULL]
dogs[, breed2 := NULL]

### 18. Add new breeds to `breeds` table ----
breeds = breeds %>%
  bind_rows(data.table(breed_name_general = (dogs %>%
                                               filter(!is.na(breed)) %>%
                                               group_by(breed) %>%
                                               summarize(n = n()) %>%
                                               filter(!breed %in% breeds$breed_name & !breed %in% breeds$breed_name_general) %>%
                                               filter(n>1) %>%
                                               pull(breed))))

### 19. Retain only valid breed strings ----
for (B in breeds$breed_name_general){
  dogs[breed == B, breed := breed]
}

### 20. Summarize breed prevalences ----
breed_prevalence = dogs %>%
  group_by(purebred,breed) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  complete(breed, purebred, fill = list(n = 0)) %>%
  group_by(breed) %>%
  mutate(breed_total = sum(n)) %>%
  ungroup() %>%
  arrange(-breed_total) %>%
  mutate(all_total = sum(n)) %>%
  mutate(breed_pct_all = breed_total/all_total)

## HEIGHT AND WEIGHT ----

### 1. Responses to Q121 and Q242 ----
dogs = dogs %>%
  merge((answers %>%
           filter(question == 121) %>%
           select(dog, height = answer) %>% 
           group_by(dog) %>% 
           summarize(height = mean(as.numeric(height), 
                                   na.rm = T))), 
        by = "dog", all.x = T) %>%
  merge((answers %>%
           filter(question == 242) %>%
           select(dog,answer) %>%
           mutate(answer = gsub("-lb"," lb", gsub("-kg", " kg", answer))) %>%
           separate(answer, into = c('value','unit'), sep = " ", convert = T) %>%
           mutate(value = abs(value)) %>%
           mutate(weight = if_else(unit == "lb",
                                   value,
                                   value*2.2)) %>%
           select(dog, weight) %>% 
           group_by(dog) %>% 
           summarize(weight = mean(as.numeric(weight), 
                                   na.rm = T))), 
        by = "dog", all.x = T) %>%
  as.data.table()

### 2. Breed-averaged heights and weights ----
breed_size = dogs %>%
  group_by(breed) %>%
  summarize(avg_height = round(mean(height, na.rm = T), digits = 0),
            avg_weight = mean(weight, na.rm = T)) %>%
  as.data.table()
breed_size[is.na(avg_height), avg_height := 2]
breed_size[is.na(avg_weight), avg_weight := (breed_size %>% filter(avg_height == 2) %>% pull(avg_weight) %>% mean(na.rm = T))]

dogs = dogs %>% merge(breed_size,
                      by = "breed",
                      all.x = T)

height_size = dogs %>%
  filter(!is.na(height)) %>%
  group_by(height) %>%
  summarize(avg_height_weight = mean(weight, na.rm = T)) %>%
  as.data.table()

dogs = dogs %>% merge(height_size,
                      by = "height",
                      all.x = T)

dogs[, height_filled := height]
dogs[, weight_filled := weight]

### 3. Fill missing heights: breed-averaged height ----
# if height (121) is NA and has breed(s), then fill mean of breed(s):
# reasoning? dog roughly size of stated breed
dogs[is.na(height) & !is.na(breed),
     height_filled := avg_height]

### 4. Fill missing heights: all-averaged height ----
# if height (121) is NA and has no breed(s), then fill 2:
# reasoning? filling average
dogs[is.na(height) & is.na(breed),
     height_filled := 2]

### 5. Fill missing weights: breed-averaged weight ----
# if weight (242) is NA and has breed(s), then fill mean of breed(s):
# reasoning? weight going to be closer to stated breed
dogs[is.na(weight) & !is.na(breed),
     weight_filled := avg_weight]

### 6. Fill missing weights: height-averaged weight ----
# if weight (242) is NA and has height (121), then fill mean of height:
# reasoning? older users may not have answered 242 but do have height
dogs[is.na(weight) & is.na(breed),
     weight_filled := avg_height_weight]

dogs = dogs %>% select(-avg_height)
dogs = dogs %>% select(-avg_weight)
dogs = dogs %>% select(-avg_height_weight)

### 7. Round off filled heights and weights ----
dogs = dogs %>%
  mutate(weight_filled = round(weight_filled, digits = 1),
         height_filled = round(height_filled, digits = 0))

## TRAINING ----

### Q#120: "Has DOG had any formal training? (check all that apply)" ----
dq.120 = answers %>%
  filter(question == 120) %>%
  select(dog,option) %>% 
  as.data.table()

dq.120[, train := 0]
dq.120[!is.na(option) & option != "Other", train := 1]
dq.120 = dq.120 %>%
  group_by(dog) %>%
  summarize(train_level = sum(train))
dogs = dogs %>%
  merge(dq.120, by = "dog", all.x = T)

## HOUSEHOLD ----

### Q#115 and Q#189: Cats in home ----
dogs = dogs %>%
  mutate(household_cats = if_else(dog %in% (
    answers %>%
      filter(question %in% c(115,189)) %>%
      select(dog,question,answer) %>%
      unnest_tokens(word, answer) %>%
      filter(word %in% c("cat",
                         "cats",
                         "feline",
                         "tabby",
                         "siamese") |
               word %like% "cat") %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE))

### Q#115 and Q#189: Dogs in household ----
dogs = dogs %>%
  mutate(household_dogs = if_else(dog %in% (
    answers %>%
      filter(question %in% c(115,189)) %>%
      select(dog,question,answer) %>%
      unnest_tokens(word, answer) %>%
      filter(word %in% c("dog","dogs","canine","mutt","terrier","lab","retriever","shepherd","collie","chihuahua","hound") | word %like% "dog") %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE))

## MEDICATIONS ----

### Q#119 "Please list any medications that DOG regularly receives." ----

# Here, we extract medication terms from Q#119 using methods in natural language processing and free text cleaning, including:
# - fixing typos and missed spaces
# - normalizing text to lowercase 
# - parsing numbers apart from letters
# - regular expressions
# - removing stop words
# - finding stem words
# - tokenizing words
# Then, matching tokens to drug names to therapeutic classifications.

# Method:
# 1. select answers
# 2. clean strings
# 3. unnest tokens
# 4. remove common and custom stop words
# 5. match token words to drug(s) via RxNorm Concept Unique Identifier (RNCUI) using search API (find 1st match)
# 6. match RNCUI to Ingredient Name (IN) and Brand Name (BN) as well as associated RNCUI for those terms
# 7. match RNCUI to WHO Anatomical Therapeutic Chemical (ATC) 1st and 2nd classifications
# 8. clean set of tokens grouped by ATC class
# 9. add any drug tokens missed by RxNorm and ATC to relevant class and save reference data for parsing again
# 10. flexibly find matches (incl. mispellings) in dog-token pairs to infer whether dog is taking drug of a given class

# https://www.nlm.nih.gov/research/umls/rxnorm/index.html
# https://www.who.int/tools/atc-ddd-toolkit/atc-classification

med_stopwords = c("vet","veterinarian",
                  "approximately",
                  "day","daily","biweekly","weekends","cycles",
                  "combination","ombination",
                  "releif","relief",
                  "scheduled","plaqu",
                  "month","monthl","monthly",
                  "dosage","dose",
                  "mg","medication","meds",
                  "0","1","2","3",
                  "supplement","holistic",
                  "repellent","spray",
                  "recurrent",
                  "heart",
                  "counter","green","comfort","laser",
                  "topical","springtime","hurts","foster",
                  "oster","chloe","question",
                  "prevantative","cooperative",
                  "finate","regime","helped",
                  "crazy","daycare","greenie","isn't",
                  "adopted","advocate")

medications = answers %>%
  filter(question == 119) %>%
  select(dog,question,answer) %>%
  mutate(answer = gsub("([a-z])([A-Z])","\\1 \\2",answer)) %>%
  unnest_tokens(word, answer) %>%
  count(word, sort = TRUE) %>%
  filter(!word %in% stop_words$word) %>%
  filter(!word %in% med_stopwords) %>%
  filter(!word %like% "mg") %>%
  filter(nchar(word) > 4) %>%
  filter(!str_detect(word, "[0-9]")) %>%
  filter(n > 1) %>%
  arrange(-n)

medications$search_rxcui = unlist(lapply(lapply(X = medications$word, FUN = function (x) {RxNormR::rx_approximateTerm(x)$approximateGroup$candidate[[1]]$rxcui}), function(x) if (length(x) == 0) 0 else x))

medications = medications %>%
  filter(search_rxcui != 0) %>%
  mutate(generic_rxcui = purrr::map_chr(.x = search_rxcui,
                                        .f = function(x) paste(rx_related_tty(x, "IN")[[1]]$conceptGroup[[1]]$conceptProperties[[1]]$rxcui, collapse = "")),
         brand_rxcui = purrr::map_chr(.x = search_rxcui,
                                      .f = function(x) paste(rx_related_tty(x, "BN")[[1]]$conceptGroup[[1]]$conceptProperties[[1]]$rxcui, collapse = "")))

medications = medications %>%
  rowwise() %>%
  mutate(therapeutic = paste(rxnorm::get_atc(generic_rxcui, "second"), collapse = "|")) %>%
  mutate(generic_name = purrr::map_chr(.x = generic_rxcui,
                                       .f = function(x) get_rx(x)))

medications = medications %>%
  filter(!is.na(therapeutic) & therapeutic != "NA" & !is.na(generic_name)) %>%
  select(word,n,rxcui=generic_rxcui,rxname=generic_name,rxatc=therapeutic)

medications = medications %>%
  separate(rxatc, into = paste("tx",1:10, sep = ""), sep = "\\|") %>%
  pivot_longer(cols = paste("tx",1:10, sep = ""),
               names_to = "txout",
               values_to = "rxatc",
               values_drop_na = T) %>%
  filter(rxatc != "NA") %>%
  select(word,n,rxcui,rxname,rxatc)

# therapeutic classes of interest:
tx_include = medications %>%
  group_by(rxatc) %>%
  summarize(n = sum(n)) %>%
  filter(n > 5) %>%
  arrange(-n) %>%
  pull(rxatc)

# drugs of interest:
px_include = medications %>%
  group_by(rxname) %>%
  summarize(n = sum(n)) %>%
  filter(n > 5) %>%
  arrange(-n) %>%
  pull(rxname)

# check how owners using words:
check_word = "monthl"
answers %>%
  filter(question==119) %>%
  filter(answer %like% check_word) %>%
  pull(answer)

# drug-tx key
med_key = answers %>%
  filter(question %in% c(119)) %>%
  select(dog,question,answer) %>%
  mutate(answer = gsub("([a-z])([A-Z])","\1 \2",answer)) %>%
  unnest_tokens(word, answer) %>%
  filter(!word %in% stop_words$word) %>%
  filter(!word %in% med_stopwords) %>%
  filter(!word %like% "mg") %>%
  filter(nchar(word) > 4) %>%
  filter(!str_detect(word, "[0-9]")) %>%
  filter(word %in% (medications %>%
                      filter(rxatc %in% tx_include | rxname %in% px_include) %>%
                      pull(word)))

med_key = medications %>% merge(med_key, by = "word")

dog_meds = answers %>%
  filter(question == 119) %>%
  select(dog,answer,notes) %>%
  merge((med_key %>%
           select(dog,word,rxcui,rxname,rxatc)), by = "dog", all = T) %>%
  arrange(dog,rxname) %>%
  filter(!is.na(rxname)) %>%
  unique() %>%
  as.data.table()

# 197 drug names
# 70 therapeutic categories

dog_meds_rxname = dog_meds %>%
  group_by(rxname) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(tot = n()) %>%
  mutate(pct = n/tot) %>%
  select(rxname,pct) %>%
  unique() %>% 
  arrange(-pct) %>%
  filter(!is.na(rxname)) %>%
  filter(rxname != "ethanol") %>%
  filter(pct > 0.005) %>%
  pull(rxname)

dog_meds_rxatc = dog_meds %>%
  group_by(rxatc) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(tot = n()) %>%
  mutate(pct = n/tot) %>%
  select(rxatc,pct) %>%
  unique() %>% 
  arrange(-pct) %>%
  filter(!is.na(rxatc)) %>%
  filter(pct > 0.0025) %>%
  pull(rxatc)

# drugs and usage in dogs
rx_dogs = c("mirtazapine",        # for appetite
            "amitriptyline",     # for behavior (general)
            "clomipramine",      # for behavior (general)
            "doxepin",           # for behavior (general)
            "fluoxetine",        # for behavior (general)
            "sertraline",        # for behavior (general)
            "paroxetine",        # for behavior (general)
            "mirtazapine",       # for behavior (general)
            "fluvoxamine",       # for behavior (general)
            "venlafaxine",       # for behavior (general)
            "bupropion",         # for behavior (general)
            "buspirone",         # for anxiety
            "diazepam",          # for anxiety
            "alprazolam",        # for anxiety
            "acepromazine",      # for anxiety
            "lorazepam",         # for anxiety
            "clorazepate",       # for anxiety
            "dexmedetomidine",   # for anxiety
            "gabapentin",        # for seizures
            "phenobarbital",     # for seizures
            "zonisamide",        # for seizures
            "levetiracetam",     # for seizures
            "clonazepam",        # for seizures
            "pregabalin",        # for seizures
            "melatonin",         # for sleep
            "benazepril",        # for heart
            "enalapril",         # for heart
            "lisinopril",        # for heart
            "tramadol",          # for pain
            "aspirin",
            "clonidine",
            "acetaminophen",
            "ivermectin",
            "praziquantel",
            "moxidectin",
            "thiabendazole",
            "fenbendazole",
            "piperazine",
            "selegiline",
            "amantadine",
            "hydrocortisone",
            "prednisone",
            "prednisolone",
            "attapulgite",
            "neomycin",
            "loratadine",
            "budesonide",
            "miconazole",
            "ondansetron",
            "scopolamine",
            "dronabinol",
            "cetirizine",
            "hydroxyzine")


drug_pain = dog_meds %>%
  filter(rxatc %in% c("antiinflammatory and antirheumatic products",
                      "analgesics")) %>%
  pull(rxname) %>% unique()

drug_parasite = dog_meds %>%
  filter(rxatc == "anthelmintics") %>%
  pull(rxname) %>% unique()

drug_heart = c("")

drug_sedative = c( "gabapentin",        # for seizures
                   "phenobarbital",     # for seizures
                   "zonisamide",        # for seizures
                   "levetiracetam",     # for seizures
                   "clonazepam",        # for seizures
                   "pregabalin")        # for seizures

drug_allergy = c("hydrocortisone",
                 "prednisone",
                 "prednisolone",
                 "attapulgite",
                 "neomycin",
                 "budesonide",
                 "miconazole",
                 "ondansetron",
                 "scopolamine",
                 "dronabinol",
                 "cetirizine",
                 "hydroxyzine",
                 "diphenhydramine",
                 "fexofenadine",
                 "loratadine",
                 "meclizine",
                 "chlorpheniramine",
                 "clemastine")

drug_psych = c("clonidine",          # off label
               "mirtazapine",        # for appetite
               "amitriptyline",     # for behavior (general)
               "clomipramine",      # for behavior (general)
               "doxepin",           # for behavior (general)
               "fluoxetine",        # for behavior (general)
               "sertraline",        # for behavior (general)
               "paroxetine",        # for behavior (general)
               "mirtazapine",       # for behavior (general)
               "fluvoxamine",       # for behavior (general)
               "venlafaxine",       # for behavior (general)
               "bupropion",         # for behavior (general)
               "buspirone",         # for anxiety
               "diazepam",          # for anxiety
               "alprazolam",        # for anxiety
               "acepromazine",      # for anxiety
               "lorazepam",         # for anxiety
               "clorazepate",       # for anxiety
               "dexmedetomidine")   # for anxiety

dogs = dogs %>%
  mutate(drug_psych = if_else(dog %in% (
    dog_meds %>%
      filter(rxname %in% drug_psych) %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE),
  drug_allergy = if_else(dog %in% (
    dog_meds %>%
      filter(rxname %in% drug_allergy) %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE),
  drug_pain = if_else(dog %in% (
    dog_meds %>%
      filter(rxname %in% drug_pain) %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE),
  drug_sedative = if_else(dog %in% (
    dog_meds %>%
      filter(rxname %in% drug_sedative) %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE))

rm(drug_allergy)
rm(drug_heart)
rm(drug_pain)
rm(drug_parasite)
rm(drug_psych)
rm(drug_sedative)

# SAVE DATA ----
## Finalize tables ----
dogs = dogs %>%
  select(dog,
         name,
         profile_public,
         image_private,
         owner,
         group,
         sex,
         neutered,
         sex_status,
         height,
         height_filled,
         weight,
         weight_filled,
         text_age,
         age,
         age_group,
         birth_date,
         birth_year,
         adopt_date,
         adopt_year,
         consent_date,
         consent_year,
         start_date,
         start_year,
         flagged_deceased_date,
         death_year,
         purebred,
         breed,
         breed.nonbreed,
         breed.alternate,
         breed1.original,
         breed2.original,
         working_dog,
         state,
         country,
         region,
         environ,
         origin,
         household_dogs,
         household_cats,
         starts_with("drug_")) %>%
  arrange(dog)

questions = questions %>% 
  filter(is.na(tags) | tags != "SPORT") %>% 
  arrange(id) %>% 
  select(id,
         string,
         tags,
         style,
         format,
         options,
         survey,
         title,
         intro)

answers = answers %>%
  select(dog,
         question,
         index = id,
         answer,
         option,
         notes,
         time,
         age,
         user) %>%
  arrange(dog,question,index)

## Save files ----
write_csv(dogs,
          paste(workDir,
                "/dat/",
                "DarwinsArk_",
                as.character(freezeDate),
                "_dogs.csv",
                sep=""))

write_csv(questions,
          paste(workDir,
                "/dat/",
                "DarwinsArk_",
                as.character(freezeDate),
                "_questions.csv",
                sep=""))

write_csv(answers,
          paste(workDir,
                "/dat/",
                "DarwinsArk_",
                as.character(freezeDate),
                "_answers.csv",
                sep=""))

## Save image ----
save.image(file = paste(workDir,
                        "/dat/", 
                        paste(paste("DarwinsArk",
                                    freezeDate, 
                                    sep="_"), 
                              ".RData", 
                              sep=""), 
                        sep=""))