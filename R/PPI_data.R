#' Building PPI network data
#'
#' @param name Set of genes named ENTREZID
#' @return PPI network data
#' @keywords Protein Protein interaction
#' @export



PPI_data = function(name){

  checkmate::assert_integer(name)

  tmp = a_human_pure[1,]
  tmp$gene_id1 = NA
  tmp$interactant_id1 = NA
  N = 1

  for (x in 1:length(name)) {
    hanghao = which(a_human_pure == name[x], arr.ind = T)

    if (dim(hanghao)[1] != 0) {
      for (y in 1:dim(hanghao)[1]) {
        tmp[N,] = a_human_pure[hanghao[y,1],]
        N = N + 1
      }
    }

  }

  # 去除只出现一次的互作

  a = tmp$gene_id1
  b = tmp$interactant_id1
  c = c(a,b)

  c = sort(table(c), decreasing = T)
  c = as.data.frame(c)

  quchu= c[c$Freq == 1,]
  quchu$c = as.integer(as.character(quchu$c))
  class(quchu$c)

  ## 备份互作数据集
  net = tmp
  hanghao = NA


  for (x in 1:length(quchu$c)) {
    hanghao = which(tmp == quchu$c[x], arr.ind = T)

    tmp = tmp[-hanghao[1],]
  }




  a = tmp$gene_id1
  b= tmp$interactant_id1
  c = unique(c(a, b))
  table(!name %in% c) # 可以看到有几个关键基因被去除了
  name[!name %in% c]  # 输出被被去除的基因[因为只有一条互作，因此也也不加了]
  # 打印到屏幕上
  if (length(name[!name %in% c]) != 0) {
    options(warn = 1)
  }
  warning(paste('These key genes have only one interaction data, so they are removed:', name[!name %in% c], '\n'))

  net = tmp
  return(net)

}

