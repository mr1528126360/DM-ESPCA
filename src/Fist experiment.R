source('DM-ESPCA.R')

data1 <- read.csv("../data/35809_28838.csv", header=TRUE)
RelationMatrix <- read.csv("../data/RelationMatrix_35809_28838.csv", header=TRUE)
Metabolic_t<- read.csv("../data/Metabolic_t.csv", header=TRUE)
Invasive_t<- read.csv("../data/Invasive_t.csv", header=TRUE)
Proliferative_t<- read.csv("../data/Proliferative_t.csv", header=TRUE)
n = 100
myString <- "start!"
print(myString)
edges = list()
gene1 = data.matrix(data1)
gene = t(gene1)
re_av = array(0,dim=c(200))
for (i in 1:1181312){
  edges[[i]] = c((RelationMatrix[i,1]+1),(RelationMatrix[i,2]+1))
}

p_t = list()
i_t = list()
m_t = list()
for (i in 1:28838)
{
  p_t[[i]] = c((Proliferative_t[i,1]))
}
for (i in 1:28838)
{
  i_t[[i]] = c((Invasive_t[i,1]))
}
for (i in 1:28838)
{
  m_t[[i]] = c((Metabolic_t[i,1]))
}
weight_list = list(m_t,p_t,i_t)

out3 = DM_ESPCA(gene, k = 3, edges, k.group=n,we=0.3, t = 0.1,niter=100, w_l=weight_list)

myString <- "end!"
print(myString)
write.csv (out3[["V"]], file ="../result/result.csv")
write.csv (out3[["U"]], file ="../result/result_u.csv")

