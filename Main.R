library(epimod)
library(epimodFBAfunctions)


wd = getwd() #Setto la working directory

setwd(wd)

#Vado a settare i valori del file da leggere (modello dell'ecoli)
model_type = "ecoli"
model_name = "iML1515"


#Prendo il riferimento al percorso della rete (GreatMode)
net_wd = paste0(wd, "/net")
net_fname = paste0(model_name, "_PN")

print(net_fname)


model = FBA4Greatmod.generation(fba_mat = './input/iML1515.mat',
                        fba_model = './net/iML1515_PN.PNPRO'
                        )

writeFBAfile(model, './result/result')

#Prove generazione + analisi 

model_gen = model.generation(net_fname = './net/iML1515_PN.PNPRO', transitions_fname = './net/transitions.cpp', fba_fname = './result/result')


model.analysis(solver_fname = 'iML1515_PN.solver', f_time = 10, s_time = 0.5)
