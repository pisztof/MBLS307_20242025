# The aim of this tutorial is to visualize protein-protein interactions observed in the study by
# Kim et al., 2021 and perform enrichment analysis on the obtained network using R. Here, we will
# perform simplified visualization compared to the original article.

# We will use the following libraries:
# Stringdb: a library to visualize Protein Protein Interactions (PPI's) based on published (experimental) data
# tidyverse: contains dplyr, tidyr and stringr: provide an easy way to filter and wrangle data
# magrittr: together with tidyverse it allows the use of %>%
# readxl and writexl: libraries to open and write microsoft excel files
# RGBL: is used in 'network' analysis using graph theory
# RGraphviz: a library to visualize PPI's
# graph and igraph: libraries that provide tools for PPI analysis

# Install and load libraries
library(STRINGdb) #to install initially BiocManager::install("STRINGdb")
string_db <- STRINGdb$new(version="11.0b", species=9606, score_threshold=200, input_directory="")
library(tidyverse)
library(readxl)
library(writexl)
library(magrittr)
library(RBGL) # BiocManager::install("RBGL")
library(Rgraphviz) # BiocManager::install("Rgraphviz")
library(graph)  # BiocManager::install("graph")
library(igraph) # from CRAN via install.packages("igraph")

# With the libraries in memory, we will start by loading in the supplementary data of the article.
# The 'raw' PPI data is in Supplementary table S3. The other tables contain the same data, however from other stages of post-processing and hence more refined.
# For the purpose of this tutorial we will use the data from Supplementary data S6 (science.abf3066_Table_S6.xlsx).
# You can find the input file in the folder "input-data" in the folder where this script was.
# This excel file contains two (2) sheets, README and "Differential Interaction Score", we need the latter.

# Load the supplementary data of S6 and the correct sheet into a new variable SuppS6
input <- read_excel("input-data/science.abf3066_Table_S6.xlsx", sheet = "Differential Interaction Score")

head(input)

# check the content of the imported data. What columns do you see? Check in README of the file what the columns mean

# Converting the data into an Graph object
# The data as it is cannot be used directly in the Graph library and its functions.
# For this, we need to turn it into a graphNEL graph.
# Basically the graph data structure are two vectors: baits and matching preys;
# They become nodes on the graph, while their interactions are displayed as edges.

# Let's first make a custom function which will allow us to create a graphNEL object.
# GraphNEL provides a very general structure for representing graphs.
# Instances of graphNEL are assumed to be simple graphs with at most one edge between any pair of nodes

makeproteingraph <- function(myBaits, myPreys)
{
  uniquenames <- unique(c(myBaits,myPreys))
  
  # Make a graph for these proteins with no edges:
  mygraph <- new("graphNEL", nodes = uniquenames)
  
  # Add edges to the graph:
  # See http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf for more examples
  mygraph2 <- addEdge(myBaits, myPreys, mygraph, 1)
  
  return(mygraph2)
}

# Create our first PPI plot.
# Now that we have a function to create graphNEL, we need the list of baits, and the list of preys
# SuppS6 contains two columns: baits and preys column.
# It also contains two additional columns, with an other identifier (uniprot) for these proteins.
# We will use baits and preys.

interactionGraph <- makeproteingraph(myBaits = input %>% select(Bait) %>% pull(),
                                     myPreys = input %>% select(Prey) %>% pull())

# call the interactionGraph variable.
# As an output you will see the number of Nodes and Edges in your graph (but not a plot yet)
interactionGraph

# to visualise the graph, we first need to define its layout and then render it, which is done with the package 'graphviz
# Read the first two pages of the vignette for the graphviz library
# https://rdrr.io/bioc/Rgraphviz/f/inst/doc/newRgraphvizInterface.pdf.

# Before you can use the renderGraph function, we need to set the layout of our graph.
# We will make use of the neato layout layoutType="neato" (there are more layout types, but 'calculating' these is slow).
# To layout your graph correctly you need the code below (just copy paste it)

interactionGraph <- layoutGraph(interactionGraph, layoutType="neato")

# now you can call rendergraph on your graphNEL-object
interactionGraph <- renderGraph(interactionGraph)

# Take a look at the graph and compare it to Figure 2A in the paper.
# It doesn't look the same but trust me it is the same PPI information as in the figure.

# Figure 2H in the article uses the same data, however now only shows the interactions
# that are commonly found in the cancer lines and not in the mcf10A line.
# In other words, they only look at those PPI's that are specific for cancer cells and not observed in the 'normal' mcf10a line.
# The supplementary data has a DIS column. The autors describe that a positive DIS means that it is specific to cancer lines (README part in TABLE S6)

sub_input <- input %>% filter(DIS > 0.5) # we limited the PPI data to ~ 369 nodes
head(sub_input)
dim(sub_input)

# make a PPI plot of the subset¶
# We have our subset, which should contain only the PPIs found in cancer lines.
# Of course we like to plot this as well, like we did with the complete dataset.

# do the following steps by yourself
# make a new graphNEL-object
# output the myGraph to see the number of nodes and edges
# layout the myGraph using the layoutGraph(.., layoutType="neato") function
# render the graph

filteredGraph <- makeproteingraph(myBaits = input %>% select(Bait) %>% pull(),
                                     myPreys = input %>% select(Prey) %>% pull()) %>% 
  layoutGraph(layoutType="neato") %>% renderGraph()


# Similar to the figure in the article, we would like to highlight the baits a bit better.
# Also the preys are a bit crowded. Below you find an example code that takes the baits,
# makes them blue, and changes the circle nodes into rectangles.
# It also changes the font size for the baits and makes the font for the preys very small (invisible).

# The details are in: https://www.bioconductor.org/packages/devel/bioc/vignettes/Rgraphviz/inst/doc/newRgraphvizInterface.pdf However for this tutorial I will provide the code (below). Try to see if you understand the code. Run it and maybe play with the code a bit to change the visualization to your liking.

# In short this is what the code below does:
# It starts by again layouting the graph (it basically gets rid of previous layouts).
# All the settings are done by creating a list.
# This list has per element a name, which is the bait or prey name,
# and per name a setting, like the fontSize, color etc.
# This allows for a per Node (which can be a bait or prey) setting.
# Here we just fill the list with one color or font size or whatever setting we want to change.
# So, all baits are blue, and fontsize 48 etc.
# For each setting we create a new list. Yes this isn't very efficient, and it can be done differently.
# But this is 'easier' to follow. This system can also be used to change the format of the EDGES.
# Settingslist (e.g. fill)
# Bait1 - lightblue
# Bait2 - lightblue
# ....
# BaitN - lightblue
# 
# SettingsList 2 (e.g font)
# Bait1 - 48
# Bait2 - 48
# ....
# BaitN - 48
# 
# etcetera
# graph.par is to set a title and subtitle
# render the graph

## https://www.bioconductor.org/packages/devel/bioc/vignettes/Rgraphviz/inst/doc/newRgraphvizInterface.pdf

sub_interactionGraph <- layoutGraph(filteredGraph, layoutType="fdp") #layoutType="fdp" is also cool, but slow, try it out!

## formatting the prey nodes (sub_input$Prey)
preyFont <- rep(1, length(sub_input$Prey)) #assign a very small font size to the prey nodes
names(preyFont) <- sub_input$Prey


## formatting the bait nodes (sub_input$Bait) 
fill <- rep("lightblue", length(sub_input$Bait)) #assign the blue color to the bait nodes
names(fill) <- sub_input$Bait

shape <- rep("box", length(sub_input$Bait)) #assign the box shape to the bait nodes
names(shape) <- sub_input$Bait

baitFont <- rep(48, length(sub_input$Bait))#assign the font size to the bait nodes
names(baitFont) <- sub_input$Bait

lwd <- rep(2, length(sub_input$Bait)) #make the borders around the bait nodes squares 2px thick.
names(lwd) <- sub_input$Bait


## the nodeRenderInfo takes these list with parameters (settings) and uses that format the graph.
nodeRenderInfo(sub_interactionGraph) <- list(fontsize=preyFont)
nodeRenderInfo(sub_interactionGraph) <- list(fontsize=baitFont) 
nodeRenderInfo(sub_interactionGraph) <- list(fill=fill)
nodeRenderInfo(sub_interactionGraph) <- list(shape=shape)
nodeRenderInfo(sub_interactionGraph) <- list(lwd=lwd)
nodeRenderInfo(sub_interactionGraph) <- list(lty=1)

## some additional markup stuff, set at title and subtitle.
graph.par(list(nodes=list(fontsize=36),
               graph=list(main="PPI NETWORK",
                          sub="DIS > 0.5",
                          cex.main=1.8, cex.sub=1.4, col.sub="gray")))

## finally we use the renderGraph function to draw the graph again
renderGraph(sub_interactionGraph)

# Maybe this plot is not 1:1 the same look and feel as Figure 2H in the article,
# but it is much more informative like this.

# It would be interesting to see if we can also learn from this graph.
# For example, one can extract information which proteins are connected to HRAS or TP53.
# For that we can use the adj(YourGraphNEL-object, "ProteinName") function. Try this function with for example TP53 as protein.

# As a sidenode: for this tutorial we used the Graph library. A more recent library is iGraph.
# You can convert Graph to iGraph and use iGraph for more advanced features.
# However, more advanced is also a bit more complex. When you have become a more advanced R user, you can switch to iGraph or the even newer ggraph library.

adj(sub_interactionGraph, "TP53")

# Another interesting thing to learn is the degree of vertices. Or in other words, how many interactions does a protein have.
# For this, we can use the function degree(GraphNEL-object) and put the output in the mydegrees variable
# Note: we have to ad graph:: before the degree function.
# This is because there is also an igraph::degree function.
# By just using degree() R doesn't know which package to use (and you get an error)

# 1. Run the graph::degree function and put the outcome in a variable: mydegrees
# This function will look at every node and determine how many connections it has and output this in a list.
# We can visualize this list in a histogram and see the distribution of nodes and connections.
# Some nodes (proteins) will have only one connection, but some will have a lot of connections.
# 2. Output the content of the mydegrees variable: use the following code, because we will sort it on the fly
# sort(mydegrees)

mydegrees <- graph::degree(sub_interactionGraph)
sort(mydegrees)
# Next we use the R hist() function to plot the histogram of this data.
hist(mydegrees)

# There is a lot of information you can extract from these networks.
# If you find this interesting, you can dive into the world of graph-theory and networks for more information.
# We will cover more advanced topics next week

# This concludes the first part of the tutorial, vizualization of the PPI networks. 
# We will now continue with using a database with known information on proteins and protein networks to extract functional annotations.


# As is often the case when using different libraries / packages, the data that we have 'as is' cannot be used by the STRING packages.
# STRING for example doesn't use the protein-protein interaction information (e.g. protein A pulled-down protein B).
# String bases protein networks on its own database of reported interactions.
# It takes your list of proteins and subsequently looks in its database if there are interactions between the other proteins in your list reported.
# Therefore we only need the list of BAITS and PREYS.
# However, String doesn't use protein names directly, it wants to use its own idendifiers, aka STRING_id.
# The data we will use is not the complete interaction data, but the data we obtained after filtering for DIS > 0.5.
# These PPIs that are specific to the cancer cells.

Baits_sup <- sub_input$Bait
Preys_sup <- sub_input$Prey
sub_interactionGraph <- makeproteingraph(Baits_sup, Preys_sup)

# All we need to do is combine Baits and Preys into one big list.
allprots <- c(Baits_sup, Preys_sup)

# Now you can proceed to put these proteins in a dataframe variable 'protein'
# and use as column name proteinID
protein <- data.frame(proteinID = allprots)

# The final step is to add a column with the STRING identifiers.
# Open the vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf
# On page 3 of the vignette (almost at the bottom) you will find an example of the MAP function.
# Use string_db$map to add the identifier column to our data.frame
# In the example in the vignette, authors create a new data.frame,
# but you can just send the output back to the protein variable : protein <- string_db$..........

protein <- string_db$map(protein, "proteinID", removeUnmappedRows = TRUE)


# don't worry if not all proteins could be mapped to a String identifier **
# we can now ask String to plot a network, based on all the proteins in our dataset.
# On page 4 you will find the string_db function to "plot a network". Try this
# If you run out of memory, reduce the size of the dataframe as in the vignette, e.g. to 500.
top_protein <- protein$STRING_id[1:500]
string_db$plot_network(top_protein)

# the obtained network is difficult to check for individual components.
# It is also very different to our PPI plots. Still these are the same protein.
# StringDB 'determines' the protein-protein interaction not based on our PPI map,
# but on what it "knows" from other experiments and based on literature.
# The colored edges between the nodes are indicative for the kind of evidence there is for such an interaction.

# As you have seen from the PPI maps, there were several clusters of proteins.
# Often such a community could be an indicator of a protein complex.
# Within string networks you can also look for such clusters.
# These clusters are now based on the text mining (literature) and observed data.
# On page 12 and 13 of the vignette, you will find an example on how String can be used to find communities (or clusters).
# Lets create a cluster list first.
# Use the function as described on page 12/13; store the ouput in the variable: clutersList
clustersList <- string_db$get_clusters(protein$STRING_id)

# Call the cluster list object to get some information on these clusters.
clustersList
length(clustersList)
lengths(clustersList)

# plot the clusters;
# if you get an error "Error in clustersList[[i]] : subscript out of bounds", fewer than X clusters were found
par(mfrow=c(4,2))
par(mar=c(1,1,1,1))
for(i in seq(1:7)){
  string_db$plot_network(clustersList[[i]])
}

# We can repeat the step above with this cluster to further take apart a cluster.
# For instance, you can cluster the first cluster another time, to go even deeper! Let;s do this.
# Use as input the clusterList[[??]] (where ?? should be replaced with the cluster number of the cluster of interest)
# and store the outcome of this clustering round in the variable: clustersList_1
# How many clusters and proteins are there (length() and lengths())
# Plot all these 'sub' clusters

clustersList_1 <- string_db$get_clusters(clustersList[[1]])
length(clustersList_1)
lengths(clustersList_1)

par(mfrow=c(3,1))
for(i in seq(1:3)){
  string_db$plot_network(clustersList_1[[i]])
}

# Now that we have seen that there are some really strong clusters present,
# we would like to know if we can determine some kind of function associated to them.
# For this, we can use Enrichment analysis (page 8), where we will use as input our initial list
# of proteins from before the clustering,
# and store the output of the analysis in the enrichment variable.
# Next, we save the content of enrichment in an excel file, which makes it easier to look at.
# Use write_xlsx(enrichment, "NameOfYourFile.xlsx")

enrichment <- string_db$get_enrichment(protein)
head(enrichment, n=20)

write_xlsx(enrichment, "enrichment-result-all.xlsx")

# Open the excel file, and see if you can derrive some cellular processes and functions that are related to cancer.
# Please note the added value of running the enrichment analysis with literature (PubMed), protein domains (pfam).

# We can do the same thing for our clusters and sub-clusters.
# As the proteins that are in a (sub) cluster are stored in the clustersList variable (and in the clustersList1 variable etc etc).
# Try it yourself.
# Generally, you will see that the enrichment analyses can find more specific functions and modules.
clusterProtein <- filter(protein,
                         STRING_id %in% clustersList)
Cluster_enrichment <- string_db$get_enrichment(clustersList[[1]])
head(Cluster_enrichment, n=20)

write_xlsx(Cluster_enrichment, "enrichment-result-cluster1.xlsx")

#### END of the tutorial ####