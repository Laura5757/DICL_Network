library(readxl)
library(tidyverse)
library(corrplot)
library(dplyr)
library(ggbeeswarm)
library(emmeans)
library(broom)
library(ggpubr)   
library(epitools)
library(GGally)
library(igraph)
library(ggnetwork)
library(intergraph)
library(xtable)
# DATA LOADING ------------------------------------------------------------
#Data reference
#https://dataverse.harvard.edu/file.xhtml?fileId=10743461&version=3.0
data <- read_csv("D:/Maria/Desktop/KTU/SA/L8/DICL_v2.csv")
View(data)
#Using lpn (Linguistic proximity index for different native languages)
df <- data[, c("country_i", "country_j","lpn", "col", "cor")] #Keep the countries, index and language sharing
#Delete redundant relationships
df <- df %>%
  rowwise() %>% 
  mutate(Node1 = min(country_i, country_j), 
         Node2 = max(country_i, country_j)) %>%
  ungroup()%>%
  distinct(Node1, Node2, .keep_all = TRUE)%>%
  select(c("country_i", "country_j","lpn", "col", "cor"))

#Filter to keep only the relationships with countries with different official or the facto languages
df <- df %>%
  filter(col == 0 & cor == 0)
#Impute zero values to NA
df$lpn[is.na(df$lpn)] <- 0 
#Filter indices greater than zero
df <- df %>%
  filter(lpn > 0)
#Filter data set for Europe and Asia
europe_asia_countries <- c(
  "Albania", "Andorra", "Belarus",
  "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Czech Republic",
  "Denmark", "Estonia", "Finland", "France", "Georgia", "Germany", "Greece",
  "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", 
  "Montenegro", "Netherlands", "North Macedonia", "Norway", "Poland", "Portugal",
  "Romania", "Russia", "Serbia", "Slovakia", "Slovenia", "Spain",
  "Sweden", "Ukraine", "United Kingdom", "India"
)
#Filtering both countries columns to keep the relationships between Asia and Europe countries languages
df_i <- df %>% 
  filter(country_i %in% europe_asia_countries)
df_j <- df_i %>%
  filter(country_j %in% europe_asia_countries)
df_j <- df_j %>%
  filter(country_i != country_j)
#3Compute the distances for later usage
df_j$distances <- 1/df_j$lpn

# DATA VISUALISATION ------------------------------------------------------
#Create network from dataframe
g_df <- graph_from_data_frame(
  d = df_j, 
  directed = FALSE 
) 
#Edges
E(g_df)$weight <- df_j$lpn
#Vertices
V(g_df)$label <- V(g_df)$name
#Plot graph
plot(g_df)

#Create filtered network for lpn greater than 0.3 
df_k <- df_j %>%
  filter(lpn > 0.3)

g <- graph_from_data_frame(
  d = df_k, 
  directed = FALSE 
) 
#Edges
E(g)$weight <- df_k$lpn
#Vertices
V(g)$label <- V(g)$name
#Plot graph
plot(g)

#Network using the distances (with filtered)
g_d <- graph_from_data_frame(
  d = df_k, 
  directed = FALSE 
) 
#Edges
E(g_d)$weight <- df_k$distances
#Vertices
V(g_d)$label <- V(g_d)$name


# CENTRALITY CALCULATION --------------------------------------------------
#Compute the degree centrality for g
degr_cent <- centr_degree(g_d, mode = 'all')
degr_cent <- degr_cent$res

#Compute the closeness centrality
clos_cent <- igraph::closeness(g_d)

#Compute betweeness centrality
betw_cent <- igraph::betweenness(g_d)

#Compute the eigenvector centrality of our network, as it uses the weigth to the conection, we use the lpn
eign_cent <- eigen_centrality(g)$vector

#Compute PageRank centrality
pr_cent <- page_rank(g)$vector



#Dataframe of the measures of centrality
measures <- data.frame(
                   degree = degr_cent, 
                   closeness = clos_cent, 
                   betweeness = betw_cent,
                   eigen = eign_cent, 
                   pagerank = pr_cent)

#By degree centrality
measures <- measures %>% arrange(desc(degree))
head(measures)
#Get the measurements table in tex format
#print(xtable(measures, type = "latex"), file = "D:/Maria/Desktop/KTU/SA/L8/measures_table.tex")

#As I am using the distances, I will use the layout Kamada Kawai because it uses longer edges to represent larger values 
plot(g_d,                          
     edge.color = 'black',
     vertex.size =degr_cent,
     vertex.shape = 'circle',      
     asp = 1,      
     layout= layout_with_kk)       

# CLIQUES CALCULATION -----------------------------------------------------
#All cliques 
length(cliques(g_d, min=1))
#Maximal cliques
count_max_cliques(g_d)
#The largest clique
largest_cliques(g_d)

#All cliques full dataframe
length(cliques(g_df, min=10))
#Maximal cliques
count_max_cliques(g_df)
largest_cliques(g_df)

# COMMUNITY DETECTION -----------------------------------------------------
## FOR THE NOT FILTERED DATAFRAME  ------------------------------------------
#Find communites using the modularity maximisation algorithm (using lpn no distances)
g_df_simplified = simplify(g_df)
df_mod_groups <- cluster_fast_greedy(g_df_simplified)
df_mod_groups <- df_mod_groups$membership
tkplot(g_df_simplified, vertex.color = df_mod_groups, 
       edge.color = 'black',
       vertex.size = centr_degree(g_df, mode = 'all')$res,         
       vertex.shape = 'circle',  
       asp = 0,
       main = 'Modularity')

#Find communites using the Girvan-Newman algorithm (using distances)
g_dd <- graph_from_data_frame(
  d = df_j, 
  directed = FALSE 
) 
E(g_dd)$weight <- df_j$distances
V(g_dd)$label <- V(g_dd)$name

dd_btw_groups <- cluster_edge_betweenness(g_dd)
dd_btw_groups <- dd_btw_groups$membership
tkplot(g_dd, vertex.color = dd_btw_groups, 
       edge.color = 'black',
       vertex.size = centr_degree(g_dd, mode = 'all')$res,         
       vertex.shape = 'circle',  
       asp = 0,
       main = 'Girvan-Newman algorithm')

#Clustering (using lpn no distances)
g_df_dendro <- cluster_fast_greedy(g_df_simplified)
plot_dendrogram(g_df_dendro)

## FOR THE FILTERED DATAFRAME----------------------------------------------
#Find communites using the modularity maximisation algorithm (using lpn no distances)
g_simplified = simplify(g)
g_mod_groups <- cluster_fast_greedy(g_simplified)
g_mod_groups <- g_mod_groups$membership
plot(g_simplified, vertex.color = g_mod_groups,
     edge.color = 'black',     
     vertex.size = centr_degree(g, mode = 'all')$res,
     vertex.shape = 'circle',  
     asp = 0,
     layout = layout_with_fr,
     main = 'Modularity')

#Find communites using the Girvan-Newman algorithm (using distances)
gd_btw_groups <- cluster_edge_betweenness(g_d) 
gd_btw_groups <- gd_btw_groups$membership
plot(g_d, vertex.color = gd_btw_groups, 
     edge.color = 'black',
     vertex.size = centr_degree(g_d, mode = 'all')$res,         
     vertex.shape = 'circle',  
     asp = 0,
     layout =  layout_with_kk,
     main = 'Girvan-Newman algorithm')

#Clustering 
g_dendro <- cluster_fast_greedy(g_simplified)
plot_dendrogram(g_dendro)




