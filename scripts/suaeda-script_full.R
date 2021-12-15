# Suaeda Landscape Genetics Analysis
# July 13, 2021
# Carina Isabella Motta

#1 LOAD PACKAGES----------------------------------------------------------------

# a vector listing package names needed for importing the DNA sequences,
#calculating genetic distance, calculated geographic distance, and performing
#a Mantel test

  package.list <- c("here", #so I don't have to deal with setting a WD
                  "adegenet", #import DNA sequences
                  "ape", #calculating genetic distance
                  "reshape2", #reshaping matrix to vector 
                  "Imap", #measuring geographic distance
                  "vegan", #mantel test
                  "tidyverse", #data cleaning
                  "dplyr", #data cleaning
                  "stringr", #data cleaning
                  "ggplot2" #mantel test visualization 
                  )

#installing the packages if they aren't already on the computer
  new.packages <- package.list[!(package.list %in% installed.packages()
                                 [,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

#and loading the packages into R with a for loop
  for(i in package.list){library(i, character.only = T)}
  
??bipartite

#2 LOADING DNA SEQUENCES--------------------------------------------------------

#exported aligned sequences from Geneious as .fasta file and edited in 
  #Notepad++ to simplify sequence names and change order 
  
#OPTION 1------------------
  #load .fasta file of aligned sequences as a DNAbin object using 'ape' package
        #read.dna(here("alignment.fasta"), format = "fasta", as.matrix = T)

#OPTION 2------------------ 
  #(use this method)
  #load .fasta file of aligned sequences as a DNAbin object using 
  #'adegenet' package
seq <- fasta2DNAbin(here("data", "2021_alignment.fasta"), quiet=FALSE, chunkSize=10, 
                     snpOnly=FALSE)

  

  #confirm it loaded correctly
  #str(seq)
  #head(seq)
  

#3 CALCULATE GENETIC DISTANCE--------------------------------------------------- 

#OPTION 1------------------
  #computes a matrix of pairwise distances from DNA sequences using different
  #models of DNA evolution  with the 'ape' package
  
  #This function resulted in a lot of 0s between samples, so I decided not 
  #to use it at first, but then I found out that the reason why it resulted
  #in more zeros than the 'dist.gene' function is because it deletes columns
  #with ambiguity codes

  ?dist.dna
  
  #I used the four models below to see if they gave different results

  gendist_raw <- 100*(dist.dna(seq, model='raw', pairwise.deletion = F, 
                         as.matrix=T))
  
  #gendist_raw_test <- 100*(dist.dna(seq_test, model='raw', pairwise.deletion = F, 
                               #as.matrix=T))

  #gendist_kimura <- 100*(dist.dna(seq, model='K80', pairwise.deletion = F, 
                               #as.matrix=T))
  
  #gendist_TN93 <- 100*(dist.dna(seq, model='TN93', pairwise.deletion = F, 
                               #as.matrix=T))
  
  #gendist_JC69 <- 100*(dist.dna(seq, model='JC69', pairwise.deletion = F, 
                               #as.matrix=T))
  
  write.csv(gendist_raw, here("data", "gendist_raw.csv"))

#OPTION 2------------------
  #computes a matrix of distances between pairs of individuals using
  #package 'ape' as well, but I do not think it deletes columns with 
  #ambiguity codes, therefore artificially inflating the number of differences
  #between samples

  ?dist.gene

  #Pairwise---
    #gendist_pairwise <- dist.gene(seq, method = "pairwise", 
                        #pairwise.deletion = FALSE,
                        #variance = FALSE)

    #write.csv(gendist_pairwise, "gendist_pairwise.csv") 
    #converted to matrix so it can be viewed 
      #gendist_pairwise <- as.matrix(gendist_pairwise)
    #confirm structure
      #str(gendist_pairwise)
    #write .csv so I could copy the sample names to make the geo_dist matrix
      #write.csv(gendist_pairwise, "gendist_pairwise.csv")
  
  #Percentage---
    #gendist_percentage <- dist.gene(seq, method = "percentage", 
                        #pairwise.deletion = FALSE,
                        #variance = FALSE)
      
#Heirfstat Package-------------
      
#the 'genet.dist' fuction from the hierfstat package is usually used in mantel 
#tests to calculate Fst values, the difference in allele frequencies on a scale
#from 0 to 1.
#however it requires a dataset of alleles from multi-locus genotypes      
      
      ?genet.dist
  

#4 MEASURE GEOGRAPHIC DISTANCES BETWEEN LOCALITIES------------------------------

  #CREATE DATAFRAME---------   
    #create dataframe of the lat/lon coordinates of each individual in the same 
    #order as the matrix of genetic distances
      #abbreviated name of locality = full locality name (number of samples)
      #ll = Las Lisas (4)
      #sca = San Carlos (3)
      #sr = Santa Rosa (4)
      #loa = Los Angeles (4)
      #laa = Las Animas (3)
      #sf = San Felipe (4)
      #si = San Ignacio (3)
      #scr = Santa Cruz (4)
      #sg = San Gregorio (1)
      
    geo_suaeda <- data.frame(name = c("1044_si", "1045_si", "1062_si", 
                                      "1064_sg", "1051_sca", "1052_sca", 
                                      "1053_sca", "1046_laa", "1047_laa", 
                                      "1048_laa", "1040_loa", "1041_loa", 
                                      "1042_loa", "1043_loa", "994_sf", 
                                      "995_sf", "1038_sf", "1039_sf", 
                                      "990_ll", "991_ll", "992_ll", 
                                      "993_ll", "1054_sr", "1055_sr", 
                                      "1056_sr", "1057_sr", "1058_scr",
                                      "1059_scr", "1060_scr", "1061_scr"
                                      ),
                          lat = c(26.8008333, 26.8008333, 26.8008333, 26.0616667,
                                  24.7997222, 24.7997222, 24.7997222, 28.8016667,
                                  28.8016667, 28.8016667, 28.9716667, 28.9716667, 
                                  28.9716667, 28.9716667,31.2880556, 31.2880556,
                                  31.2880556, 31.2880556, 31.7841667, 31.7841667,
                                  31.7841667, 31.7841667, 28.9666667, 28.9666667,
                                  28.9666667, 28.9666667, 28.8, 28.8, 28.8, 28.8),
                          lon = c(-113.1480556, -113.1480556, -113.1480556,
                                  -112.272777777777, -112.1152778, -112.1152778,
                                  -112.1152778, -113.3313889, -113.3313889,
                                  -113.3313889, -113.5466667, -113.5466667,
                                  -113.5466667, -113.5466667, -114.8927778,
                                  -114.8927778, -114.8927778, -114.8927778,
                                  -114.6288889, -114.6288889, -114.6288889,
                                  -114.6288889, -112.1416667, -112.1416667, 
                                  -112.1416667, -112.1416667, -111.9166667,
                                  -111.9166667,-111.9166667,-111.9166667))
      
  #CREATE FUNCTIONS TO MAKE MATRIX---------     
    #the code below was borrowed from a Data Mining PowerPoint Presentation
    #that can be found at http://www2.cs.uh.edu/~ceick/UDM/UDM/Rfunc.pdf 
    #(slides 14 - 16) using the gdist function from the 'Imap' package and 
    #creates a function that generates a matrix of geographic distances
      
      ?gdist

      ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
        
        if(nrow(m) !=ncol(m)) stop("Supplied matrix must be square.")
        if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
        else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
        m[tri] <- t(m)[tri]
        return(m)
      }

  GeoDistanceInMetersMatrix <- function(df.geopoints){
  #returns a matrix (M) of distances between geographic points
  #M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i])
  #and (df.geopoints$lat[j], df.geopoints$lon[j])
  #The row and column names are given by df.geopoints$name
  
    GeoDistanceInMeters <- function(g1, g2){
    #Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    #The 1st value in the returned vector is the distance between g1[[1]]
    #and g2[[1]]. The 2nd value in the returned vector is the distance between 
    #g1[[2]] and g2[[2]]. Etc. Each g1[[x]] or g2[[x]] must be a list with the 
    #named elements "index", "lat", and "lon". 
    #E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), 
    #list("index"=3, "lat" = 12.1, "lon" = 13.2))
    
      DistM <- function(g1, g2){
        require("Imap")
        return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, 
                                                    lon.1=g1$lon,
                                                    lat.2=g2$lat, 
                                                    lon.2=g2$lon,
                                                    units="m")))
    
       }
      return(mapply(DistM, g1, g2))
      
      }
  
    n.geopoints <- nrow(df.geopoints)
  
    df.geopoints$index <- 1:n.geopoints
  
    list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, 
                       function(x){return(list(x))})
  
    mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, 
                                                      list.geopoints,
                                                      GeoDistanceInMeters), 
                                                      "lower")
  
    rownames(mat.distances) <- df.geopoints$name
    colnames(mat.distances) <- df.geopoints$name
  
    return(mat.distances)
  }

  #CREATE MATRIX
    geo_dist <- round(GeoDistanceInMetersMatrix(geo_suaeda) / 1000)
    
    #divide by 1000 to make units km


#5 CONDUCT MANTEL TEST----------------------------------------------------------

#using package 'vegan'
    
    ?mantel

#Mantel r values can fall within a range between -1 to 1. 
#An r value of -1 suggests a strong negative correlation, 
#0 suggests no relationship at all and 
#1 suggests a strong positive relationship.

#I ran mantel tests using the "pearson" method with each genetic distance 
#matrix I calculated to see how the models I used differed
    
mantel(gendist_raw, geo_dist, method ="pearson", permutations = 9999)
    #Mantel statistic r: 0.2025 
    #Significance: 0.0175 

#mantel(gendist_raw, geo_dist, method="pearson", permutations = 9999)
    #Mantel statistic r: 0.201 
    #Significance: 0.0207

#mantel(gendist_TN93, geo_dist, method="pearson", permutations = 9999)
    #Mantel statistic r: 0.2065 
    #Significance: 0.0164

#mantel(gendist_JC69, geo_dist, method="pearson", permutations = 9999)
    #Mantel statistic r: 0.2025 
    #Significance: 0.0201

#All models resulted in a Mantel statistic r value 0.201 - 0.2065 
#(weak positive correlation) and all being significant with p < 0.05

#Because the different models for calculating genetic distance resulted in 
#more or less the same mantel test results, I decided to use the "raw" model
#going forward

#COASTAL DISTANCE

  coast <- read.csv(here("data", "geodist_coastline.csv"), row.names = 1, header = T,
                    stringsAsFactors = F)
  
    coast[ , 1:30] <- apply(coast[ , 1:30], 2,            
                        function(x) as.numeric(as.character(x)))
    
    geodist_coastline <- coast %>% rename_all(~stringr::str_replace(.,"^X",""))
    
    geodist_coast <- as.matrix(geodist_coastline)
  
  mantel(gendist_raw, geodist_coast, method = "pearson", permutations = 9999)

#SAME COAST
  
  same <- read.csv(here("data", "geodist_samecoast.csv"), row.names = 1, header = T,
                    stringsAsFactors = F)
  
  same[ , 1:30] <- apply(same[ , 1:30], 2,            
                          function(x) as.numeric(as.character(x)))
  
  geodist_samecoast <- same %>% rename_all(~stringr::str_replace(.,"^X",""))
  
  geodist_same <- as.matrix(geodist_samecoast)
  
  mantel(gendist_raw, geodist_same, method = "pearson", permutations = 9999)

#OCEAN DISTANCE
  
  ocean <- read.csv(here("data","geodist_oceandist.csv"), row.names = 1, header = T,
                    stringsAsFactors = F)
  
  ocean[ , 1:30] <- apply(ocean[ , 1:30], 2,            
                          function(x) as.numeric(as.character(x)))
  
  geodist_oceandist <- ocean %>% rename_all(~stringr::str_replace(.,"^X",""))
  
  geodist_ocean <- as.matrix(geodist_oceandist)
  
  mantel(gendist_raw, geodist_ocean, method = "pearson", permutations = 9999)


#BARRIER
  
  barrier <- read.csv(here("data","geodist_barrier.csv"), row.names = 1, header = T,
                    stringsAsFactors = F)
  
  #made contents of data.frame characters for boxplot
  barrier[ , 1:30] <- apply(barrier[ , 1:30], 2,            
                          function(x) as.character(as.character(x)))
  
  geodist_barrier <- barrier %>% rename_all(~stringr::str_replace(.,"^X",""))
  
  geodist_barrier <- as.matrix(geodist_barrier)
  
  str(geodist_barrier)
  mantel(gendist_raw, geodist_barrier, method = "pearson", permutations = 9999)


#6 VISUALIZATION---------------------------------------------------------------

?as.vector

#convert matrices to vectors
gend <- as.vector(gendist_raw)
geod <- as.vector(geo_dist)
geoc <- as.vector(geodist_coast)
geoo <- as.vector(geodist_ocean)
geob <- as.vector(geodist_barrier)

#make new data frame with vectorized distance matrices
mat <- data.frame(gend, geod)
mat2 <- data.frame(gend, geoc)
mat3 <- data.frame(gend, geoo)
mat4 <- data.frame(gend, geob)


#plot euclidean geographic distance
vis <- ggplot(mat, aes(y = gend, x = geod)) + 
  geom_smooth(method = "lm", alpha =0.51, colour = "blue") + 
  geom_point(size =2, alpha = 1) +
  labs(y = "Genetic Distance (%)", x = "Geographic Distance (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 10), 
         axis.text.y = element_text(face = "bold", size = 10, colour = "black"), 
         axis.title= element_text(face = "bold", size = 12, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         strip.text = element_text(size = 9, face = "bold"))

vis

#plot coastal distance
vis2 <- ggplot(mat2, aes(y = gend, x = geoc)) + 
  geom_smooth(method = "lm", alpha =0.51, colour = "blue") + 
  geom_point(size =2, alpha = 1) +
  labs(y = "Genetic Distance (%)", x = "Coastal Distance (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 10), 
         axis.text.y = element_text(face = "bold", size = 10, colour = "black"), 
         axis.title= element_text(face = "bold", size = 12, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         strip.text = element_text(size = 9, face = "bold"))

vis2

#plot oceanic distance 
vis3 <- ggplot(mat3, aes(y = gend, x = geoo)) + 
  geom_smooth(method = "lm", alpha =0.51, colour = "blue") + 
  geom_point(size =2, alpha = 1) +
  labs(y = "Genetic Distance (%)", x = "Oceanic Distance (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 10), 
         axis.text.y = element_text(face = "bold", size = 10, colour = "black"), 
         axis.title= element_text(face = "bold", size = 12, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         strip.text = element_text(size = 9, face = "bold"))

vis3

#plot with barrier of peninsula absent/present 
vis4 <- ggplot(mat4, aes(y = gend, x = geob)) + 
  geom_boxplot() +
  labs(y = "Genetic Distance (%)", x = "Biogeographic Barrier Absent/Present (0/1)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 10), 
         axis.text.y = element_text(face = "bold", size = 10, colour = "black"), 
         axis.title= element_text(face = "bold", size = 12, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         strip.text = element_text(size = 9, face = "bold"))
vis4

#plot with barrier of peninsula absent/present with points
vis5 <- ggplot(mat4, aes(y = gend, x = geob)) + 
  geom_boxplot() +
  geom_jitter(alpha = 0.5, width = 0.2) +
  labs(y = "Genetic Distance (%)", x = "Biogeographic Barrier Absent/Present (0/1)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 10), 
         axis.text.y = element_text(face = "bold", size = 10, colour = "black"), 
         axis.title= element_text(face = "bold", size = 12, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         strip.text = element_text(size = 9, face = "bold"))
vis5


#export plots as pdf
ggsave(
  filename = here::here("final-figures", "euclidean.pdf"),
  plot = vis,
  device = "pdf",
  scale = 1,
  dpi = 720,
  limitsize = TRUE,
)

ggsave(
  filename = here::here("final-figures", "coastal.pdf"),
  plot = vis2,
  device = "pdf",
  scale = 1,
  dpi = 720,
  limitsize = TRUE,
)


ggsave(
  filename = here::here("final-figures", "oceanic.pdf"),
  plot = vis3,
  device = "pdf",
  scale = 1,
  dpi = 720,
  limitsize = TRUE,
)


ggsave(
  filename = here::here("final-figures", "barrier.pdf"),
  plot = vis4,
  device = "pdf",
  scale = 1,
  dpi = 720,
  limitsize = TRUE,
)

ggsave(
  filename = here::here("final-figures", "barrier_points.pdf"),
  plot = vis5,
  device = "pdf",
  scale = 1,
  dpi = 720,
  limitsize = TRUE,
)
