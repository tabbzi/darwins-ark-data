scores = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_factor-scores.csv")
  
dogs = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_dogs.csv")

anc = read_csv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno/dat/adm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_global-adm_supervised_K-115.csv")

scores.dogs = scores %>%
  merge(dogs,
        by = "dog",
        all.x = T)

scores.summ = scores.dogs %>%
  select(factor,breed,score_raw) %>%
  group_by(breed,factor) %>%
  summarize(score.mean = mean(score_raw, na.rm = T),
            score.sd = sd(score_raw, na.rm = T))

scores.all.summ = scores.dogs %>%
  select(factor,breed,score_raw) %>%
  group_by(factor) %>%
  summarize(score.mean = mean(score_raw, na.rm = T),
            score.sd = sd(score_raw, na.rm = T))

scores.all.norms = scores.all.summ %>%
  group_by(factor) %>%
  summarize(val = rnorm(10000,
                        mean = score.mean,
                        sd = score.sd)) %>%
  ungroup() %>%
  mutate(breed = "all dogs")

scores.norms = scores.summ %>%
  filter(breed %in% (anc %>% filter(pct > 0.85) %>% group_by(pop) %>% summarize(n = n()) %>% filter(n>10) %>% pull(pop))) %>%
  group_by(breed,factor) %>%
  summarize(val = rnorm(10000,
                        mean = score.mean,
                        sd = score.sd)) %>%
  ungroup() %>%
  group_by(factor,breed)

# tests
ks.test.list = list()
i=1
for (fa in unique(scores.norms$factor)){
  for (B in unique(scores.norms$breed)){
    this.test = ks.test(x= (scores.norms %>%
                                      filter(factor == fa) %>%
                                      filter(breed == B) %>%
                                      pull(val)),
                                y= (scores.all.norms %>%
                                      filter(factor == fa) %>%
                                      pull(val)))
    
    ks.test.list[[i]] = data.table(factor=fa,
                                   breed=B,
                                   p=this.test$p.value)
    
    i=i+1
  }
}
ks.test.list = rbindlist(ks.test.list)

ks.test.list= ks.test.list %>%
  group_by(factor) %>%
  mutate(p.adj = p.adjust(p = p, method = "BH"))

#scores.norms = bind_rows(scores.norms,
#                         scores.all.norms)

scores.norms = scores.norms %>%
  merge(ks.test.list, by = c("factor","breed"))

scores.norms = scores.norms %>%
  mutate(sig = p.adj < 0.05)

p = ggplot() +
  geom_density(data = scores.all.norms,
               aes(x = val),
               fill = "#505050",
               alpha = 0.50) +
  geom_density(data = scores.norms,
               aes(x = val,
                   group = breed,
                   color = p.adj < 0.05),
               alpha=0.25) +
  facet_wrap(.~factor,
             scales = "free_x",
             ncol = 5,
             nrow = 5) +
  theme_pubr() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


ggsave(filename = paste(workDir,
                        "/gwa/phe/",
                        "fa.breed-agg.sample-dists.pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 3 *3,
       height = 2 * 3)



####
scores.dogs = scores %>%
  merge(dogs,
        by = "dog",
        all.x = T) %>%
  group_by(factor) %>%
  mutate(val = scale(score_raw)) %>%
  ungroup()

scores.dogs.suff = scores.dogs %>%
  filter(breed %in% (scores.dogs %>%
                       filter(!is.na(breed)) %>%
                       group_by(breed,factor) %>%
                       summarize(n=n()) %>%
                       filter(n>=100) %>%
                       pull(breed)))

ks.test.list = list()
i=1
for (fa in unique(scores.dogs.suff$factor)){
  for (B in unique(scores.dogs.suff$breed)){
    this.test = ks.test(x= (scores.dogs.suff %>%
                              filter(factor == fa) %>%
                              filter(breed == B) %>%
                              pull(val)),
                        y= (scores.dogs %>%
                              filter(factor == fa) %>%
                              pull(val)))
    
    ks.test.list[[i]] = data.table(factor=fa,
                                   breed=B,
                                   p=this.test$p.value)
    
    i=i+1
  }
}
ks.test.list = rbindlist(ks.test.list)

ks.test.list= ks.test.list %>%
  mutate(p.adj = p.adjust(p = p, method = "BH"))

scores.dogs = scores.dogs %>%
  merge(ks.test.list, by = c("factor","breed"))

p = ggplot() +
  geom_density(data = scores.dogs,
               aes(x = val),
               fill = "#505050",
               alpha = 0.50) +
  geom_density(data = scores.dogs,
               aes(x = val,
                   group = breed,
                   color = p.adj < 0.05),
               alpha=0.25) +
  facet_wrap(.~factor,
             scales = "free",
             ncol = 5,
             nrow = 5) +
  scale_color_discrete(type = rev(greenfocus[1:2])) +
  theme_pubr()


ggsave(filename = paste("~/Dropbox (UMass Medical School)/Thesis/figures/",
                        "fa.by-breed.dists.pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 26,
       height = 14)



dogs %>% 
  group_by(purebred) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n))


# dim reduct scores
install.packages('uwot')
library(uwot)

fa.mat = scores %>% select(dog,factor,score_fill_mean) %>%
  pivot_wider(id_cols = dog,
              names_from = factor,
              values_from = score_fill_mean) %>%
  select(-dog) %>%
  as.matrix()

fa.umap = umap(na.omit(fa.mat))
colnames(fa.mat)

fa.plot = bind_cols(as.data.frame(fa.umap),
                    as.data.frame(na.omit(fa.mat)))


fa.plot %>%
  pivot_longer(cols = -c(V1,V2),
               names_to = "factor",
               values_to = "score") %>%
  group_by(factor) %>%
  mutate(score = scale(score)) %>%
  arrange(abs(score)) %>%
  ungroup() %>%
  mutate(factor = as.integer(factor)) %>%
  ggplot(aes(x = V1,
             y = V2,
             color = score)) +
  facet_wrap(.~factor) +
  geom_point(alpha = 0.5,
             size = 0.2) +
  scale_color_viridis_c() +
  theme_pubclean()
