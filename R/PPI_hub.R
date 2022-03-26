#' Using algorithm to calculate Hub gene
#'
#' @param PPI_data_symbol A set of interaction data
#' @return The result of hub gene found by the algorithm
#' @keywords Protein Protein interaction Hub gene
#' @importFrom igraph graph_from_data_frame as_ids V diameter vcount distances neighbors components gsize degree betweenness closeness evcent eccentricity
#' @importFrom utils head
#' @export



PPI_hub = function(PPI_data_symbol){
  # 得到互作网络
  net = igraph::graph_from_data_frame(PPI_data_symbol, directed = FALSE)

  node = igraph::as_ids(igraph::V(net))   # 将节点转化为字符串的格式
  diam = igraph::diameter(net)    # 计算图的直径
  result = node
  result =as.data.frame(result)   # result用来存放算法结果的临时数据

  ## Radiality 径向度算法 ******************************************************************************************
  {
    n = igraph::vcount(net)
    result$Radiality = NA

    for (y in 1:length(node)) {
      distance =igraph::distances(net, node[y], to = V(net))
      sum = 0

      for (x in 1:length(distance)) {

        if (distance[x] != 0) {
          a = diam + 1 - distance[x]
          sum = sum + a
        }


      }

      Rad = sum / (n - 1)
      result[y,2] = Rad

    }


    Rad = head(result[order(result[,2], decreasing = T), ], n = 10)
    data_Radiality = Rad$result


  }


  ## MNC 最大邻居连通度度算法 ******************************************************************************************
  {
    result$MNC = NA
    result$DMNC = NA
    n = length(node)

    for (x in 1:n) {
      neighbor_node = igraph::neighbors(net, node[x])    ## 指定节点的邻居节点
      neighbor_node = as_ids(neighbor_node)
      # 第2小步：将当前节点A和他的邻居节点构成限定子图
      # 将不在邻居节点中的节点帅选出来
      node_no = node   ## 构建邻居中不存在的节点原始数据框

      for (y in 1:n) {
        node_no[y] = !(node[y] %in% neighbor_node)
      }


      node_no = as.logical(node_no)
      node_no_exit = node[node_no]   ## 里面存放的是邻居限定子图中的节点
      ## 第二步构建邻居限定子图

      net_tmp = net - node_no_exit
      #  plot(net_tmp)

      ## 第三步：找出最大连通分支，并计算其中的点数
      igraph::components(net_tmp)

      result[x,3] = max(igraph::components(net_tmp)$csize)

      side = gsize(net_tmp)
      result[x,4] = side / max(igraph::components(net_tmp)$csize)

    }

    MNC = head(result[order(result[,3], decreasing = T),], n = 10)
    data_MNC = MNC$result
  }

  ## DMNC 最大邻居连通密度算法 ******************************************************************************************
  {
    DMNC = head(result[order(result[,4], decreasing = T),], n = 10)
    data_DMNC = DMNC$result
  }


  ## ClusteringCoefficient 集聚系数算法 ******************************************************************************************
  {
    n = igraph::vcount(net)  # 获取节点数
    N = n * (n-1)
    N   # 这些节点之间可能存在的边数
    result$CCT = NA
    for (x in 1:n) {
      road = length(neighbors(net, node[x]))
      CCT = 2 * road / N
      result[x,5] = CCT
    }


    CCT = head(result[order(result[,5], decreasing = T), ], n = 10)
    data_ClusteringCoefficient = CCT$result
  }


  # degree 点度中心度算法 ********************************************************************************************
  {
    igraph::V(net)$deg = igraph::degree(net, mode = 'all')
    data_deg = head(igraph::V(net)[order(igraph::V(net)$deg, decreasing = T)], n = 10)
    data_deg = as_ids(data_deg)
  }

  ## Betweenness 中间中心度算法  ***************************************************************************************
  {
    igraph::V(net)$bte = igraph::betweenness(net, directed = F)
    data_bte =  head(igraph::V(net)[order(igraph::V(net)$bte, decreasing = T)], n = 10)
    data_bte = as_ids(data_bte)
  }


  ## closeness 接近中心度算法 ******************************************************************************************
  {

    igraph::V(net)$clo = igraph::closeness(net, mode = 'in')    # in为无向图

    data_clo = head(igraph::V(net)[order(igraph::V(net)$clo, decreasing = T)], n = 10)
    data_clo = as_ids(data_clo)
  }

  ## Eigenvector centrality 特征向量中心度算法 ******************************************************************************************
  {
    igraph::V(net)$EC = igraph::evcent(net, scale = F)$vector
    data_EC = head(igraph::V(net)[order(igraph::V(net)$EC, decreasing = T)], n = 10)
    data_EC = as_ids(data_EC)
  }

  ## Eccentricity算法 偏心率中心度算法 ******************************************************************************************
  {
    igraph::V(net)$Eccentricity = igraph::eccentricity(net)
    data_Eccentricity = head(igraph::V(net)[order(igraph::V(net)$Eccentricity, decreasing = T)], n = 10)
    data_Eccentricity = as_ids(data_Eccentricity)
  }




  data = cbind(data_Radiality, data_MNC, data_DMNC, data_ClusteringCoefficient,
               data_deg, data_bte, data_clo, data_EC, data_Eccentricity)
  return(data)
}





