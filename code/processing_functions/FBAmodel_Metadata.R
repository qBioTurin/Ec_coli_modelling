
FBAmodel.metadata <- function(model, model.name, wd, bigg.path) {
  
  br = ExtractEx(model, model.name = model.name, bigg.path = bigg.path)
  
  data = rbind(data.frame(value = br$lb,
                          dir = rep("lowbnd", length(br$lb)),
                          react.id = br$react.id,
                          react.index = br$index),
               data.frame(value = br$ub,
                          dir = rep("uppbnd", length(br$ub)),
                          react.id = br$react.id,
                          react.index = br$index))
  
  all_react = data.frame(React_ID = model@react_id,
                         ReactionType = rep(".", dim(model@S)[2]),
                         React_Lb = model@lowbnd,
                         React_Ub = model@uppbnd,
                         ReactionDistrict = rep(".", dim(model@S)[2]),
                         ReactionPos = 1:dim(model@S)[2])
  
  all_react$ReactionDistrict[data$react.index] = "boundary"
  all_react$ReactionDistrict[which(all_react$ReactionDistrict != "boundary")] = "core"
  
  a = dplyr::filter(all_react, grepl("EX_", React_ID))
  a$ReactionType = rep("Exchange", length(a$ReactionType))
  
  b = dplyr::filter(all_react, grepl("DM_", React_ID))
  b$ReactionType = rep("Demand/Sink", length(b$ReactionType))
  
  c = dplyr::filter(all_react, grepl("sink_", React_ID))
  c$ReactionType = rep("Demand/Sink", length(c$ReactionType))
  
  d = dplyr::filter(all_react, grepl("trans", React_ID))
  d$ReactionType = rep("Transcription", length(d$ReactionType))
  
  e = dplyr::filter(all_react, grepl("repl", React_ID))
  e$ReactionType = rep("Replication", length(e$ReactionType))
  
  f = dplyr::filter(all_react, grepl("pbios", React_ID))
  f$ReactionType = rep("Biosynthesis", length(f$ReactionType))
  
  g = dplyr::filter(all_react, !grepl("BIOMASS|biomass|rep|tra|EX_|DM_|sink_|pbios", React_ID))
  g$ReactionType = rep("Internals/Transporters", length(g$ReactionType))
  
  h = all_react[which(model@obj_coef == 1), ]
  h$ReactionType = rep("Objective", length(h$ReactionType))
  
  i = dplyr::filter(all_react, grepl("BIOMASS", React_ID))
  i = i[which(i$React_ID != h$React_ID) ,]
  i$ReactionType = rep("Secondary_biomass", length(i$ReactionType))
  
  all_react = rbind(a, b, c, d, e, f, g, h, i)
  
  all_react$React_ID = gsub("\\(", replacement = "_", all_react$React_ID)
  all_react$React_ID = gsub("\\[", replacement = "_", all_react$React_ID)
  all_react$React_ID = gsub("\\)", replacement = "", all_react$React_ID)
  all_react$React_ID = gsub("\\]", replacement = "", all_react$React_ID)
  
  all_react$React_ID = gsub("_c_", replacement = "_c", all_react$React_ID)
  all_react$React_ID = gsub("__", replacement = "_", all_react$React_ID)
  
  return(list(all_react, br))
  
}
