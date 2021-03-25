InitialMotilitySelect <- function(ntrack, domain, domain_x, Xmax) {
  
  Index <- matrix(1:prod(dim(domain)),dim(domain)) #index position of ith site
  Index[domain != 1] <- NA #identify red cells
  
  Index <- Index[!is.na(Index)]
  ValidX <- domain_x[Index] > 0 & domain_x[Index] < Xmax
  Index <- cbind(Index[ValidX], domain_x[Index[ValidX]], domain_x[Index[ValidX]] < Xmax/2)
  IndexLHS <- Index[order(Index[, 2], decreasing = TRUE),]
  IndexLHS <- IndexLHS[IndexLHS[,3] == 1,]
  IndexRHS <- Index[order(Index[, 2], decreasing = FALSE),]
  IndexRHS <- IndexRHS[IndexRHS[,3] == 0,]
  
  CellSelectedIndex <- c(IndexLHS[1:ceiling(ntrack/2),1], IndexRHS[1:ceiling(ntrack/2),1])[1:ntrack]
  RowPosCell <- ((CellSelectedIndex - 1) %% nrow(domain)) + 1
  ColPosCell <- ceiling(CellSelectedIndex/nrow(domain))
  
  data <- as.data.frame(cbind(CellSelectedIndex, RowPosCell, ColPosCell))
  
  return(data)
  
}
