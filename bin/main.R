#### PREAMBLE ####

# SET DIRECTORY
workDir = "."
setwd(workDir)

# LOAD LIBRARIES
require(tidyverse)
require(data.table)
require(zoo)
require(zipcodeR)
require(anytime)
require(lubridate)
require(stringi)
require(stringr)
require(stringdist)
require(english)
require(tm)
require(topicmodels)
require(FactoMineR)
require(lattice)
require(tidytext)
require(rxnorm) # devtools::install_github("nt-williams/rxnorm")
require(RxNormR) # devtools::install_github("mpancia/RxNormR")
data(stop_words)

# SET DATES
freezeDate = "YYYYMMDD"
freezeDir = paste(workDir,"raw",freezeDate, sep = "/")

#### LOAD DATA ####
# LOAD table `dogs`:
dogs_raw = read_csv(paste(freezeDir,"dogs.csv", sep = "/"),
                    col_names = TRUE,
                    na = c("NA", "NULL", "", "null")) %>% as.data.table()

# select columns from `dogs`:
dogs = dogs_raw[fake == 0,
                list(id,
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
                     owner_id = owner,
                     group_id = group_id,
                     profile_public = is_public,
                     image_id = image,
                     image_private = image_private,
                     working_dog = is_working
                )]

# LOAD table `formats`:
formats = read_csv(paste(freezeDir,"formats.csv", sep = "/"),
                   col_names = TRUE,
                   na = c("NA", "NULL", "", "null")) %>% as.data.table()

# prepare survey options from `formats`:
formats.options = list()
i = 1
for (j in formats$id){
  dt = data.table(id = j,
                  index = as.character(seq(0,length(unlist(strsplit(formats[id == j]$options, "|", fixed = T)))-1)),
                  option = unlist(strsplit(formats[id == j]$options, "|", fixed = T)))
  formats.options[[i]] = dt
  i = i + 1
}
formats.options = rbindlist(formats.options) %>% 
  as.data.table() %>% 
  na.omit()

# load table `surveys`:
surveys = read_csv(paste(freezeDir,"surveys.csv", sep = "/"),
                   col_names = TRUE,
                   na = c("NA", "NULL", "", "null")) %>% as.data.table()

# load table `questions`:
questions = read_csv(paste(freezeDir,"questions.csv", sep = "/"),
                     col_names = TRUE,
                     na = c("NA", "NULL", "", "null")) %>%
  merge(formats, by.x = "format", by.y = "id", all.x = T) %>%
  merge(surveys, by.x = "survey", by.y = "id", all.x = T) %>%
  as.data.table()

# load table `answers`:
answers_raw = read_csv(file = paste(freezeDir,"answers.csv", sep = "/"),
                       col_names = TRUE,
                       na = c("NA", "NULL", "null"))

# parse multi-option rows from `answers`:
answers = answers_raw %>%
  as.data.table() %>%
  mutate(answer = if_else(answer == "",
                          NA_character_,
                          answer)) %>%
  filter(!is.na(answer)) %>%
  select(-id) %>%
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

answers[!question %in% (questions %>% filter(format %in% c(3,4)) %>% pull(id)), index := as.character(answer)]
answers[question %in% c(111,112), index := as.character(nest)]

# set age values for rows in `answers`:
answers = answers %>%
  merge(formats.options,
        by.x = c("format","index"),
        by.y = c("id","index"),
        all.x = T) %>%
  mutate(date = anydate(time)) %>%
  merge((dogs %>% select(id, birth_date)), by.x = "dog", by.y = "id") %>%
  mutate(age_weeks = difftime(date, birth_date, units = "weeks")) %>%
  mutate(age_weeks = replace(age_weeks, age_weeks<=0, NA)) %>%
  mutate(age_years = time_length(interval(birth_date,date), "years")) %>%
  mutate(age_years = replace(age_years, age_years<=0, NA)) %>%
  arrange(dog,question,user,time)

# LOAD table `zipcodes`:
zips = read_csv(file = paste(freezeDir,"zipcodes.csv", sep = "/"),
                na = c("NA", "NULL", "", "null")) %>% as.data.table()
colnames(zips) = c("owner","zipcode")

# PARSE zipcodes:
zips[, zipcode := gsub("\\-.*", "", zipcode)] # get only zipcodes preceding '-'
zips[nchar(zipcode) == 4, zipcode := str_pad(zipcode, width=5, side="left", pad="0")] # if 4 digits, add 0 in front of zipcode
zips[nchar(zipcode) != 5, zipcode := NA] # if still not 5 digits, remove zipcode

# GET zipcode ~ state data:
zips = merge(zips,zip_code_db, by = "zipcode", all.x = T)

# DEFINE regions:
zips[state == 'AK', region := "Noncontiguous (Alaska)"]
zips[state == 'HI', region := "Noncontiguous (Hawaii)"]
zips[state == 'PR', region := "Noncontiguous (Puerto Rico)"]
zips[state == 'AS', region := "Noncontiguous (American Samoa)"]

zips[state %in% c("CT", "ME", "MA", "NH", "RI", "VT"), region := "Northeast (New England)"]
zips[state %in% c("NJ", "NY", "PA"), region := "Northeast (Mid-Atlantic)"]
zips[state %in% c("DE", "FL", "GA", "MD", "NC", "SC", "VT", "DC", "WV", "VA"), region := "South (South Atlantic)"]
zips[state %in% c("AL", "KY", "MS", "TN"), region := "South (East South Central)"]
zips[state %in% c("AR", "LA", "OK", "TX"), region := "South (West South Central)"]
zips[state %in% c("IL", "IN", "MI", "OH", "WI"), region := "Midwest (East North Central)"]
zips[state %in% c("IA", "KS", "MN", "MO", "NE", "ND", "SD"), region := "Midwest (West North Central)"]
zips[state %in% c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY"), region := "West (Mountain)"]
zips[state %in% c("CA", "OR", "WA"), region := "West (Pacific)"]

# GET census data:
census = read_csv(file = paste(workDir,"ref","US_census_2010.csv", sep="/")) %>% as.data.table()
colnames(census) = c("zipcode", "population", "landmass", "density")

# DEFINE densities:
urban = subset(census, density >= 700)
suburban = subset(census, density >= 100 & density < 700)
rural = subset(census, density < 100)

# DEFINE environ based on densities:
zips[zipcode %in% urban$zipcode, environ := "urban"]
zips[zipcode %in% suburban$zipcode, environ := "suburban"]
zips[zipcode %in% rural$zipcode, environ := "rural"]

# ASSIGN zipcode,state,region,environ to `dogs`:
dogs = dogs %>%
  merge((zips %>% select(owner,zipcode,state,region,environ)),
        by.x = "owner_id",
        by.y = "owner",
        all.x = T) %>% 
  arrange(id) %>%
  as.data.table()

# FILL default `environ` as 'suburban':
dogs[, environ_filled := environ]
dogs[is.na(environ), environ_filled := "suburban"]

# ASSIGN origin to `dogs`:
dogs = merge(dogs,
             answers[question == 117, .SD[which.max(time)], by = "dog"][,c("dog","option")],
             by.x = "id",
             by.y = "dog",
             all.x = T) %>% setnames("option","origin")

# FILL default `origin` as 'Other':
dogs[, origin_filled := origin]
dogs[is.na(origin), origin_filled := "Other"]

#### DATES AND TIMES ####
# REPLACE WITH EARLIEST ANSWER TIMESTAMP:
start_dates = subset((
  answers_raw %>%
    group_by(dog) %>%
    arrange(-desc(time)) %>%
    slice(1)
), select = c("dog","time")) %>% as.data.table()

colnames(start_dates) = c("dog","start_date")
start_dates[, start_date := anydate(start_date)]

dogs = merge(dogs, start_dates, all.x = T, by.x = "id", by.y = "dog")
setkey(dogs, id)

#### BREEDS AND REGISTERED PUREBRED STATUS ####
# table `breeds`
breeds = read_csv(file = paste(freezeDir,"breeds.csv", sep = "/"),
                  na = c("NA", "NULL", "null")) %>% as.data.table()
breeds[, breed_name := tolower(breed_print)]

breeds.list = c(breeds$id)
names(breeds.list) = breeds$breed_name

# 1 = registered purebred, 0 = not registered, -1 = unknown
dogs[, purebred := tolower(purebred)]
dogs[purebred == "1", purebred := "yes"]
dogs[purebred != "yes" & purebred != "no", purebred := NA]

# Synchronize Breed Identifiers (NOTE: best attempts made to synchronize old free-entry breeds to new options, not error-free)
# save original entries
dogs[, breed1_historical := breed1]
dogs[, breed2_historical := breed2]

# if encoding problems:
dogs[, breed1 := tolower(iconv(enc2utf8(breed1), sub = "byte"))]
dogs[, breed2 := tolower(iconv(enc2utf8(breed2), sub = "byte"))]

dogs[breed1 %in% c("","?","null","NULL"), breed1 := NA]
dogs[breed2 %in% c("","?","null","NULL"), breed2 := NA]

dogs[breed1 == "german sherpherd dog", breed1 := "german shepherd dog"]
dogs[breed2 == "german sherpherd dog", breed2 := "german shepherd dog"]

dogs[breed1 == "german shepherd", breed1 := "german shepherd dog"]
dogs[breed2 == "german shepherd", breed2 := "german shepherd dog"]

dogs[breed1 == "pit", breed1 := "pitbull"]
dogs[breed2 == "pit", breed2 := "pitbull"]

dogs[breed1 == "lab", breed1 := "labrador"]
dogs[breed2 == "lab", breed2 := "labrador"]

dogs[like(breed1,"mix"), purebred := "no"]
dogs[like(breed1,"mix"), breed1 := NA]

dogs[like(breed2,"mix"), purebred := "no"]
dogs[like(breed2,"mix"), breed2 := NA]

dogs[like(breed1,"american eskimo"), breed1 := "american eskimo dog"]
dogs[like(breed2,"american eskimo"), breed2 := "american eskimo dog"]

dogs[like(breed1,"bassett"), breed1 := "basset hound"]
dogs[like(breed2,"bassett"), breed2 := "basset hound"]
dogs[like(breed1,"bassett hound"), breed1 := "basset hound"]
dogs[like(breed2,"bassett hound"), breed2 := "basset hound"]


dogs[like(breed1,"berger blanc suisse"), breed1 := "berger blanc suisse shepherd"]
dogs[like(breed2,"berger blanc suisse"), breed2 := "berger blanc suisse shepherd"]

dogs[like(breed1,"biewer"), breed1 := "biewer terrier"]
dogs[like(breed2,"biewer"), breed2 := "biewer terrier"]

dogs[like(breed1,"black & tan coon hound"), breed1 := "black and tan coonhound "]
dogs[like(breed2,"black & tan coon hound"), breed2 := "black and tan coonhound "]
dogs[like(breed1,"black & tan coonhound"), breed1 := "black and tan coonhound "]
dogs[like(breed2,"black & tan coonhound"), breed2 := "black and tan coonhound "]

dogs[like(breed1,"black mouth curr"), breed1 := "black mouth cur"]
dogs[like(breed2,"black mouth curr"), breed2 := "black mouth cur"]

dogs[like(breed1,"heeler"), breed1 := "australian cattle dog"]
dogs[like(breed2,"heeler"), breed2 := "australian cattle dog"]
dogs[like(breed1,"cattle dog"), breed1 := "australian cattle dog"]
dogs[like(breed2,"cattle dog"), breed2 := "australian cattle dog"]

dogs[like(breed1,"bouvier"), breed1 := "bouvier des flandres"]
dogs[like(breed2,"bouvier"), breed2 := "bouvier des flandres"]

dogs[like(breed1,"canaan"), breed1 := "canaan dog"]
dogs[like(breed2,"canaan"), breed2 := "canaan dog"]

dogs[like(breed1,"dutch shepherd sheepdog"), breed1 := "dutch shepherd"]
dogs[like(breed2,"dutch shepherd sheepdog"), breed2 := "dutch shepherd"]

dogs[like(breed1,"fox red lab"), breed1 := "labrador retriever"]
dogs[like(breed2,"fox red lab"), breed2 := "labrador retriever"]

dogs[like(breed1,"german wire hair"), breed1 := "german wirehaired pointer"]
dogs[like(breed2,"german wire hair"), breed2 := "german wirehaired pointer"]

dogs[like(breed1,"king charles cavalier"), breed1 := "cavalier king charles spaniel"]
dogs[like(breed2,"king charles cavalier"), breed2 := "cavalier king charles spaniel"]
dogs[like(breed1,"king charles spaniel"), breed1 := "cavalier king charles spaniel"]
dogs[like(breed2,"king charles spaniel"), breed2 := "cavalier king charles spaniel"]

dogs[like(breed1,"lagotti romagnoli"), breed1 := "lagotto romagnolo"]
dogs[like(breed2,"lagotti romagnoli"), breed2 := "lagotto romagnolo"]
dogs[like(breed1,"lagotto romagnolo"), breed1 := "lagotto romagnolo"]
dogs[like(breed2,"lagotto romagnolo"), breed2 := "lagotto romagnolo"]
dogs[like(breed1,"lagotto romangnolo"), breed1 := "lagotto romagnolo"]
dogs[like(breed2,"lagotto romangnolo"), breed2 := "lagotto romagnolo"]

dogs[like(breed1,"llewellin"), breed1 := "llewellin setter"]
dogs[like(breed2,"llewellin"), breed2 := "llewellin setter"]
dogs[like(breed1,"llewellyn setter"), breed1 := "llewellin setter"]
dogs[like(breed2,"llewellyn setter"), breed2 := "llewellin setter"]
dogs[like(breed1,"llewelyn setter"), breed1 := "llewellin setter"]
dogs[like(breed2,"llewelyn setter"), breed2 := "llewellin setter"]

dogs[like(breed1,"mountain curr"), breed1 := "mountain cur"]
dogs[like(breed2,"mountain curr"), breed2 := "mountain cur"]

dogs[like(breed1,"nz hunterway"), breed1 := "huntaway"]
dogs[like(breed2,"nz hunterway"), breed2 := "huntaway"]

dogs[like(breed1,"patterdale"), breed1 := "patterdale terrier"]
dogs[like(breed2,"patterdale"), breed2 := "patterdale terrier"]

dogs[like(breed1,"red tick coonhound"), breed1 := "redtick coonhound"]
dogs[like(breed2,"red tick coonhound"), breed2 := "redtick coonhound"]

dogs[like(breed1,"spinone"), breed1 := "spinoni italiani"]
dogs[like(breed2,"spinone"), breed2 := "spinoni italiani"]
dogs[like(breed1,"spinone italiano"), breed1 := "spinoni italiani"]
dogs[like(breed2,"spinone italiano"), breed2 := "spinoni italiani"]

dogs[like(breed1,"taiwanese (formosan) mountain dog"), breed1 := "taiwan dog"]
dogs[like(breed2,"taiwanese (formosan) mountain dog"), breed2 := "taiwan dog"]
dogs[like(breed1,"thaiwaneese mountian dog"), breed1 := "taiwan dog"]
dogs[like(breed2,"thaiwaneese mountian dog"), breed2 := "taiwan dog"]

dogs[like(breed1,"wire hair fox terrier"), breed1 := "wire fox terrier"]
dogs[like(breed2,"wire hair fox terrier"), breed2 := "wire fox terrier"]
dogs[like(breed1,"wire haired fox terrier"), breed1 := "wire fox terrier"]
dogs[like(breed2,"wire haired fox terrier"), breed2 := "wire fox terrier"]
dogs[like(breed1,"wirehaired fox terrier"), breed1 := "wire fox terrier"]
dogs[like(breed2,"wirehaired fox terrier"), breed2 := "wire fox terrier"]

dogs[like(breed1,"working cocker"), breed1 := "cocker spaniel"]
dogs[like(breed2,"working cocker"), breed2 := "cocker spaniel"]
dogs[like(breed1,"working type cocker spaniel"), breed1 := "cocker spaniel"]
dogs[like(breed2,"working type cocker spaniel"), breed2 := "cocker spaniel"]

# nonbreeds: village, indigenous, and landrace dogs
# pariah, village, indigenous

dogs[, nonbreed := F]

dogs[like(breed1,"africanis"), nonbreed := T]
dogs[like(breed2,"africanis"), nonbreed := T]

dogs[like(breed1,"potcake"), nonbreed := T]
dogs[like(breed2,"potcake"), nonbreed := T]
dogs[like(breed1,"pot cake"), nonbreed := T]
dogs[like(breed2,"pot cake"), nonbreed := T]

dogs[like(breed1,"village"), nonbreed := T]
dogs[like(breed2,"village"), nonbreed := T]

dogs[like(breed1,"pariah"), nonbreed := T]
dogs[like(breed2,"pariah"), nonbreed := T]

dogs[like(breed1,"sato"), nonbreed := T]
dogs[like(breed2,"sato"), nonbreed := T]

dogs[like(breed1,"local"), nonbreed := T]
dogs[like(breed2,"local"), nonbreed := T]

dogs[like(breed1,"street"), nonbreed := T]
dogs[like(breed2,"street"), nonbreed := T]


# if 'breed1' or 'breed2' exactly matches a 'breed_name' in 'breeds'
# then set breed1_proper = T
dogs[,breed1_proper := F]
dogs[,breed2_proper := F]

for (B in breeds$breed_name){
  dogs[breed1 == B, breed1_proper := T]
  dogs[breed2 == B, breed2_proper := T]
}

# if 'breed1_proper; != T and 'breed1' or 'breed2' %like$ a 'breed_name' in 'breeds'
# then set breed1 = B

for (B in breeds$breed_name){
  dogs[breed1_proper == F & like(breed1, B), breed1 := B] 
  dogs[breed2_proper == F & like(breed2, B), breed2 := B]  
}

dogs[breed1_id == 0, breed1_id := NA]
dogs[breed2_id == 0, breed2_id := NA]

# For NEW Dogs without `breed1` and `breed2`, match _id to breeds
# dogs[is.na(breed1) & !is.na(breed1_id), breed1 := breeds[id == breed1_id, breed_name]]

breed1_id_match = merge(dogs[is.na(breed1) & !is.na(breed1_id)], breeds, by.x = "breed1_id", by.y = "id", all.x = T)[, c("id","breed1_id","breed_name")]
setkey(breed1_id_match,id)
setkey(dogs,id)
dogs[is.na(breed1) & !is.na(breed1_id)]$breed1 <- breed1_id_match$breed_name

breed2_id_match = merge(dogs[is.na(breed2) & !is.na(breed2_id)], breeds, by.x = "breed2_id", by.y = "id", all.x = T)[, c("id","breed2_id","breed_name")]
setkey(breed2_id_match,id)
setkey(dogs,id)
dogs[is.na(breed2) & !is.na(breed2_id)]$breed2 <- breed2_id_match$breed_name

# For OLD Dogs without `breed1_id` and `breed2_id`, match (now corrected) breeds to _id
breed1_id_match = merge(dogs[!is.na(breed1) & is.na(breed1_id)], breeds, by.x = "breed1", by.y = "breed_name", all.x = T)[, c("id.x","id.y","breed1")][, .(id = id.x, breed1_id = id.y, breed1 = breed1)]
setkey(breed1_id_match,id)
setkey(dogs,id)
dogs[!is.na(breed1) & is.na(breed1_id)]$breed1_id <- breed1_id_match$breed1_id

breed2_id_match = merge(dogs[!is.na(breed2) & is.na(breed2_id)], breeds, by.x = "breed2", by.y = "breed_name", all.x = T)[, c("id.x","id.y","breed2")][, .(id = id.x, breed2_id = id.y, breed2 = breed2)]
setkey(breed2_id_match,id)
setkey(dogs,id)
dogs[!is.na(breed2) & is.na(breed2_id)]$breed2_id <- breed2_id_match$breed2_id

# if `breed1` is NA but `breed2` is not, switch
replace_dogs = dogs[is.na(breed1) & !is.na(breed2)]$id
dogs[id %in% replace_dogs, breed1 := breed2]
dogs[id %in% replace_dogs, breed1_id := breed2_id]
dogs[id %in% replace_dogs, breed2 := NA]
dogs[id %in% replace_dogs, breed2_id := NA]

# final conversions of breed1 / breed2 text (for dogs with breed ids)
# Elinor's conversions:
dogs[breed1 %like% "golden" & breed1 %like% "doodle", breed1 := "goldendoodle"]
dogs[breed2 %like% "golden" & breed2 %like% "doodle", breed2 := "goldendoodle"]
dogs[breed1 %like% "labra" & breed1 %like% "doodle", breed1 := "labradoodle"]
dogs[breed2 %like% "labra" & breed2 %like% "doodle", breed2 := "labradoodle"]


dogs[breed1 %like% "dachshund" & breed1 != "dachshund", alt.breed1 := breed1]
dogs[breed2 %like% "dachshund" & breed2 != "dachshund", alt.breed2 := breed2]
dogs[breed1 %like% "dachshund", breed1 := "dachshund"]
dogs[breed2 %like% "dachshund", breed2 := "dachshund"]

dogs[breed1 %like% "chihuahua" & breed1 != "chihuahua", alt.breed1 := breed1]
dogs[breed2 %like% "chihuahua" & breed2 != "chihuahua", alt.breed2 := breed2]
dogs[breed1 %like% "chihuahua", breed1 := "chihuahua"]
dogs[breed2 %like% "chihuahua", breed2 := "chihuahua"]

dogs[breed1 %like% "standard poodle", alt.breed1 := "standard poodle"]
dogs[breed2 %like% "standard poodle", alt.breed2 := "standard poodle"]
dogs[breed1 == "standard poodle", breed1 := "poodle"]
dogs[breed2 == "standard poodle", breed2 := "poodle"]

dogs[breed1 %like% "miniature poodle", alt.breed1 := "miniature poodle"]
dogs[breed2 %like% "miniature poodle", alt.breed2 := "miniature poodle"]
dogs[breed1 == "miniature poodle", breed1 := "toy poodle"]
dogs[breed2 == "miniature poodle", breed2 := "toy poodle"]

# conversions for owner-reported breed accuracy
dogs[breed1 == "rough collie", alt.breed1 := "rough collie"]
dogs[breed2 == "rough collie", alt.breed2 := "rough collie"]
dogs[breed1 == "rough collie", breed1 := "collie"]
dogs[breed2 == "rough collie", breed2 := "collie"]

dogs[breed1 == "smooth collie", alt.breed1 := "smooth collie"]
dogs[breed2 == "smooth collie", alt.breed2 := "smooth collie"]
dogs[breed1 == "smooth collie", breed1 := "collie"]
dogs[breed2 == "smooth collie", breed2 := "collie"]

dogs[breed1 == "miniature australian shepherd", alt.breed1 := "miniature australian shepherd"]
dogs[breed2 == "miniature australian shepherd", alt.breed2 := "miniature australian shepherd"]
dogs[breed1 == "miniature australian shepherd", breed1 := "australian shepherd"]
dogs[breed2 == "miniature australian shepherd", breed2 := "australian shepherd"]

dogs[breed1 == "parson russell terrier", alt.breed1 := "parson russell terrier"]
dogs[breed2 == "parson russell terrier", alt.breed2 := "parson russell terrier"]
dogs[breed1 == "parson russell terrier", breed1 := "jack russell terrier"]
dogs[breed2 == "parson russell terrier", breed2 := "jack russell terrier"]

dogs[breed1 == "russell terrier", alt.breed1 := "russell terrier"]
dogs[breed2 == "russell terrier", alt.breed2 := "russell terrier"]
dogs[breed1 == "russell terrier", breed1 := "jack russell terrier"]
dogs[breed2 == "russell terrier", breed2 := "jack russell terrier"]

dogs[breed1 == "flat-coated retrievers", breed1 := "flat-coated retriever"]
dogs[breed2 == "flat-coated retrievers", breed2 := "flat-coated retriever"]

dogs = dogs %>% select(-c("breed1_proper","breed2_proper"))

#### SEX AND SPAY/NEUTER ####
dogs[, sex := tolower(sex)]

dogs[, neutered := tolower(neutered)]
dogs[neutered != "yes" & neutered != "no", neutered := NA]

dogs[sex == "male" & neutered == "yes", sex_status :=  "neutered male"]
dogs[sex == "male" & neutered == "no", sex_status :=  "intact male"]
dogs[sex == "male" & is.na(neutered), sex_status :=  "unspecified male"]
dogs[sex == "female" & neutered == "yes", sex_status :=  "spayed female"]
dogs[sex == "female" & neutered == "no", sex_status :=  "intact female"]
dogs[sex == "female" & is.na(neutered), sex_status :=  "unspecified female"]

#### BODY SIZE ####
dogs = dogs %>%
  merge((answers %>%
           filter(question == 121) %>%
           select(dog, height = answer) %>% 
           group_by(dog) %>% 
           summarize(height = mean(as.numeric(height), 
                                   na.rm = T))), 
        by.x = "id", by.y = "dog", all.x = T) %>%
  merge((answers %>%
           filter(question == 242) %>%
           select(dog,answer) %>%
           mutate(answer = gsub("-lb"," lb", gsub("-kg", " kg", answer))) %>%
           separate(answer, into = c('value','unit'), sep = " ", convert = T) %>%
           mutate(value = abs(value)) %>%
           mutate(weight = if_else(unit == "lb",
                                   as.numeric(value),
                                   as.numeric(value)*2.2)) %>%
           select(dog, weight) %>% 
           group_by(dog) %>% 
           summarize(weight = mean(as.numeric(weight), 
                                   na.rm = T))), 
        by.x = "id", by.y = "dog", all.x = T) %>%
  as.data.table()

# breed means:
breed_size = dogs %>%
  group_by(breed1) %>%
  summarize(avg_height = round(mean(height, na.rm = T), digits = 0),
            avg_weight = mean(weight, na.rm = T)) %>%
  as.data.table()
breed_size[is.na(avg_height), avg_height := 2]
breed_size[is.na(avg_weight), avg_weight := (breed_size %>% filter(avg_height == 2) %>% pull(avg_weight) %>% mean(na.rm = T))]

dogs = dogs %>% merge(breed_size, by = "breed1", all.x = T)

height_size = dogs %>%
  filter(!is.na(height)) %>%
  group_by(height) %>%
  summarize(avg_height_weight = mean(weight, na.rm = T)) %>%
  as.data.table()

dogs = dogs %>% merge(height_size, by = "height", all.x = T)

dogs[, height_filled := height]
dogs[, weight_filled := weight]

# if height (121) is NA and has breed(s), then fill mean of breed(s):
# reasoning? dog must be around the size of stated breed
dogs[is.na(height) & !is.na(breed1), height_filled := avg_height]

# if height (121) is NA and has no breed(s), then fill 2:
# reasoning? filling average
dogs[is.na(height) & is.na(breed1), height_filled := 2]

# if weight (242) is NA and has breed(s), then fill mean of breed(s):
# reasoning? weight going to be closer to stated breed
dogs[is.na(weight) & !is.na(breed1), weight_filled := avg_weight]

# if weight (242) is NA and has height (121), then fill mean of height:
# reasoning? older users may not have answered 242 but do have height
dogs[is.na(weight) & is.na(breed1), weight_filled := avg_height_weight]

dogs = dogs %>% select(-avg_height)
dogs = dogs %>% select(-avg_weight)
dogs = dogs %>% select(-avg_height_weight)

setkey(dogs,id)

dogs = dogs %>%
  select(dog=id,
         name,
         profile_public,
         image_private,
         owner=owner_id,
         sex,
         neutered,
         sex_status,
         height,
         height_filled,
         weight,
         weight_filled,
         birth_date,
         consent_date,
         start_date,
         flagged_deceased_date,
         purebred,
         nonbreed,
         breed1,
         breed1_id,
         breed1_historical,
         alt.breed1,
         breed2,
         breed2_id,
         breed2_historical,
         alt.breed2,
         working_dog,
         zipcode,
         state,
         region,
         environ,
         environ_filled,
         origin,
         origin_filled)

#### HOUSEHOLD DATA ####
# Here, goal is to get a rough, normalized size of household and presence or absence of other animals (cats, dogs) from free text response, not a precise count of household members.

nppl = answers %>%
  filter(question == 114) %>%
  select(dog,answer) %>%
  mutate(result = answer) %>%
  mutate(result = str_replace(result, ",\ and\ ", ",")) %>%
  mutate(result = str_replace(result, "\ and\ ", ",")) %>%
  mutate(result = str_replace(result, "&", ",")) %>%
  mutate(result = str_replace(result, ",\ ", ",")) %>%
  mutate(result = gsub("[^0-9,]", "", result)) %>%
  as.data.table()

# parse into rough household size scaled from 0 to 4 (Likert)
nppl[str_length(result) == "", parsed := 1.0]
nppl[str_length(result) == 2, parsed := 1.0]
nppl[str_length(result) == 3, parsed := 1.5] # might be 1 or 2
nppl[str_length(result) == 4, parsed := 2.0]
nppl[str_length(result) == 5, parsed := 2.5] # might be 2 or 3
nppl[str_length(result) == 6, parsed := 3.0]
nppl[str_length(result) == 7, parsed := 3.5] # might be 3 or 4
nppl[str_length(result) == 8, parsed := 4] 
nppl[str_length(result) == 9, parsed := 4.5] # might be 4 or 5
nppl[str_length(result) == 10, parsed := 5]
nppl[str_length(result) == 11, parsed := 5.5]
nppl[str_length(result) == 12, parsed := 6]
nppl[str_length(result) == 13, parsed := 6.5]
nppl[str_length(result) == 14, parsed := 7]
nppl[str_length(result) == 15, parsed := 7.5]
get_nppl = function (x) {
  as.numeric(length(str_split(x, ",", simplify = T)[str_split(x, ",", simplify = T) != ""]))
}
nppl = nppl %>%
  rowwise() %>%
  mutate(parsed = if_else(result %like% ",",
                          get_nppl(result),
                          parsed)) %>%
  as.data.table()
nppl[is.na(parsed), parsed := 1.5]

nppl = nppl %>%
  merge((dogs %>% select(dog)), by = "dog", all.y = T) %>%
  mutate(parsed = if_else(is.na(parsed),
                          1.5,
                          parsed)) %>%
  mutate(household_scale_ppl = round(scales::rescale(scale(parsed, center = T, scale = T), to = c(0,4)))) %>%
  select(dog,household_scale_ppl)

dogs = dogs %>%
  merge(nppl, by = "dog", all.x = T) %>%
  mutate(household_scale_ppl = if_else(is.na(household_scale_ppl),
                                       2,
                                       household_scale_ppl))

# pets in household

# household_cats
dogs = dogs %>%
  mutate(household_cats = if_else(dog %in% (
    answers %>%
      filter(question %in% c(115,189)) %>%
      select(dog,question,answer) %>%
      unnest_tokens(word, answer) %>%
      filter(word %in% c("cat","cats","feline","tabby")) %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE))

# household_dogs
dogs = dogs %>%
  mutate(household_dogs = if_else(dog %in% (
    answers %>%
      filter(question %in% c(115,189)) %>%
      select(dog,question,answer) %>%
      unnest_tokens(word, answer) %>%
      filter(word %in% c("dog","dogs","canine","terrier","lab","retriever","shepherd","collie","chihuahua")) %>%
      pull(dog) %>%
      unique()
  ),
  TRUE,
  FALSE))

# Q120, formal training
training = answers %>%
  filter(question == 120) %>%
  select(dog,option) %>% as.data.table()
training[, train := 0]
training[!is.na(option) & option != "Other", train := 1]
training = training %>%
  group_by(dog) %>%
  summarize(train_level = sum(train))
dogs = dogs %>%
  merge(training, by = "dog", all.x = T)

#### MEDICATION DATA ####
# Q#119 "Please list any medications that DOG regularly receives."

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
  filter(question %in% c(119)) %>%
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

test_atc = "antiepileptics"

dog_meds %>%
  filter(rxatc == test_atc) %>%
  pull(rxname) %>%
  unique()

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

#### ALLERGENS ####
# Q#183 ""

# Here, we extract medication terms from Q#183 using methods in natural language processing and free text cleaning, including:
# - fixing typos and missed spaces
# - normalizing text to lowercase 
# - parsing numbers apart from letters
# - regular expressions
# - removing stop words
# - finding stem words
# - tokenizing words
# Then, matching tokens to food allergens.

#https://appliednetsci.springeropen.com/articles/10.1007/s41109-018-0109-9

library(tm)
library(stringi)
library(proxy)
library(corpus)
library(FactoMineR)
require(stringr)
require(stringdist)
require(english)
require(tm)
require(topicmodels)
require(lattice)
require(tidytext)
data(stop_words)

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
allergen_pca_dog = allergen_pca$ind$coord %>% as.data.table()
allergen_pca_dog$dog = answers %>%
  filter(question %in% c(183)) %>%
  pull(dog)

allergen_pca_dog = allergen_pca_dog %>%
  merge((answers %>% filter(question == 181) %>% select(dog,itch=answer)), by = "dog", all.x = T) %>%
  merge((answers %>% filter(question == 182) %>% select(dog,vomit=answer)), by = "dog", all.x = T) %>%
  merge((answers %>% filter(question == 178) %>% select(dog,sus_allergy=answer)), by = "dog", all.x = T)

#### AGE ####
dogs = dogs %>% 
  merge(
    (answers %>%
       group_by(dog) %>%
       summarize(age_min = min(age_years, na.rm = T),
                 age_max = max(age_years, na.rm = T),
                 age_mean = mean(age_years, na.rm = T),
                 age_group = if_else(mean(age_years, na.rm = T) >= 1.5,
                                     "adult",
                                     "puppy",
                                     "adult")) %>%
       filter(!is.nan(age_mean))),
    all.x = T,
    by = "dog") %>%
  as.data.table()

#### SAVE DATA #### 
save.image(file = paste(workDir, "dat", paste(paste("DarwinsArk", freezeDate, sep="_"), ".RData", sep=""), sep="/"))

write_csv((dogs %>% select(dog,
                           name,
                           owner,
                           sex,
                           neutered,
                           sex_status,
                           birth_date,
                           consent_date,
                           start_date,
                           flagged_deceased_date,
                           age_min,
                           age_max,
                           age_mean,
                           age_group,
                           purebred,
                           breed1_parsed = breed1,
                           breed2_parsed = breed2,
                           breed1_inputted = breed1_historical, 
                           breed2_inputted = breed2_historical,
                           breed1_alternate = alt.breed1,
                           breed2_alternate = alt.breed2,
                           nonbreed,
                           working_dog,
                           region,
                           environ,
                           origin,
                           height,
                           weight,
                           height_filled,
                           weight_filled,
                           drug_psych,
                           drug_sedative,
                           drug_pain,
                           drug_allergy) %>% arrange(dog)),
          paste(workDir, "dat", paste(paste("DarwinsArk", freezeDate, "dogs", sep="_"), ".csv", sep=""), sep="/"))

write_csv((answers %>% 
             arrange(dog,question)), 
          paste(workDir, "dat", paste(paste("DarwinsArk", freezeDate, "answers", sep="_"), ".csv", sep=""), sep="/"))

write_csv((questions %>% 
             filter(is.na(tags) | tags != "SPORT") %>% 
             arrange(id) %>% 
             select(id,string,tags,style,format,options,survey,title,intro)), 
          paste(workDir, "dat", paste(paste("DarwinsArk", freezeDate, "questions", sep="_"), ".csv", sep=""), sep="/"))
