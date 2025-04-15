# ----------------------------
# Enhanced GRN Visualization with Split MoR_sum and Clean Layout
# ----------------------------

# ---- Load Required Packages ----
if (!require("igraph")) install.packages("igraph", dependencies = TRUE)
if (!require("visNetwork")) install.packages("visNetwork", dependencies = TRUE)
if (!require("RColorBrewer")) install.packages("RColorBrewer", dependencies = TRUE)

library(igraph)
library(visNetwork)
library(RColorBrewer)

# ---- USER INPUT ----
user_seed <- 100
num_nodes <- 30
num_edges <- 60

# ---- Create Random Directed Graph ----
set.seed(user_seed)
g <- make_empty_graph(n = num_nodes, directed = TRUE)
V(g)$name <- as.character(1:num_nodes)

edge_list <- data.frame(from = character(), to = character(), type = character())
while (nrow(edge_list) < num_edges) {
  from <- sample(V(g)$name, 1)
  to <- sample(setdiff(V(g)$name, from), 1)
  if (!are.connected(g, from, to)) {
    type <- sample(c("upregulation", "inhibition"), 1)
    edge_list <- rbind(edge_list, data.frame(from = from, to = to, type = type))
  }
}
g <- add_edges(g, t(as.matrix(edge_list[, c("from", "to")])))
E(g)$type <- edge_list$type

# ---- Community Detection ----
communities <- cluster_infomap(g)
V(g)$community <- communities$membership

# ---- Centrality and Degree ----
btw <- betweenness(g, directed = TRUE)
in_deg <- degree(g, mode = "in")
out_deg <- degree(g, mode = "out")

# ---- MoR_sum Calculation (Split In/Out) ----
mor_in <- rep(0, vcount(g))
mor_out <- rep(0, vcount(g))

for (i in seq_along(V(g))) {
  node <- V(g)[i]
  
  # Incoming
  incoming_edges <- incident(g, node, mode = "in")
  for (e in incoming_edges) {
    if (E(g)[e]$type == "upregulation") {
      mor_in[i] <- mor_in[i] + 1
    } else if (E(g)[e]$type == "inhibition") {
      mor_in[i] <- mor_in[i] - 1
    }
  }
  
  # Outgoing
  outgoing_edges <- incident(g, node, mode = "out")
  for (e in outgoing_edges) {
    if (E(g)[e]$type == "upregulation") {
      mor_out[i] <- mor_out[i] + 1
    } else if (E(g)[e]$type == "inhibition") {
      mor_out[i] <- mor_out[i] - 1
    }
  }
}

mor_score <- mor_in + mor_out  # still used for node sizing

# ---- Assign Active/Inactive ----
set.seed(user_seed)
activity_status <- sample(c(1, -1), vcount(g), replace = TRUE)
V(g)$status <- activity_status

# ---- Highlight Max Centrality ----
V(g)$highlight <- FALSE
for (comm in unique(V(g)$community)) {
  comm_nodes <- which(V(g)$community == comm)
  max_node <- comm_nodes[which.max(btw[comm_nodes])]
  V(g)$highlight[max_node] <- TRUE
}

# ---- Color Setup ----
n_comms <- length(unique(V(g)$community))

if (!require("viridis")) install.packages("viridis")
library(viridis)
comm_palette <- viridis(n_comms, option = "K")


hex_to_rgba <- function(hex, alpha = 1) {
  rgb <- col2rgb(hex) / 255
  sprintf("rgba(%d,%d,%d,%.2f)", rgb[1]*255, rgb[2]*255, rgb[3]*255, alpha)
}

node_colors <- sapply(1:vcount(g), function(i) {
  base_col <- comm_palette[V(g)$community[i]]
  alpha <- 1  # fully opaque
  hex_to_rgba(base_col, alpha)
})

status_label <- ifelse(V(g)$status == 1, "Active", "Inactive")
status_border <- ifelse(V(g)$status == 1, "blue", "red")

# ---- Nodes ----
nodes <- data.frame(
  id = V(g)$name,
  label = paste0("Gene ", V(g)$name,
                 "\nIn: ", in_deg,
                 " | Out: ", out_deg,
                 " | Btwn_cen: ", round(btw, 2),
                 "\nMoR_in: ", mor_in,
                 " | MoR_out: ", mor_out,
                 "\nStatus: ", status_label),
  group = V(g)$community,
  color.background = node_colors,
  color.border = status_border,
  color.highlight.background = "#FFD700",
  color.highlight.border = "black",
  font = list(size = 20, face = "Arial", color = "black", bold = TRUE),
  shape = ifelse(V(g)$highlight, "star", "dot"),
  value = mor_score,
  shadow = TRUE
)

# ---- Edges ----
edges <- as_data_frame(g, what = "edges")
edges <- data.frame(
  from = edges$from,
  to = edges$to,
  arrows = "to",
  color = ifelse(edges$type == "inhibition", "red", "blue"),
  dashes = ifelse(edges$type == "inhibition", TRUE, FALSE),
  smooth = TRUE
)
# 2. Apply arrowhead shape customization
edges$arrows <- lapply(edges$color, function(clr) {
  if (clr == "red") {
    list(to = list(enabled = TRUE, type = "bar"))
  } else {
    list(to = list(enabled = TRUE, type = "arrow"))
  }
})
# ---- Interactive Plot ----
gg <- visNetwork(nodes, edges, height = "700px", width = "100%", background = "white") %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    nodesIdSelection = list(enabled = TRUE, useLabels = TRUE)
  ) %>%
  visIgraphLayout(layout = "layout_with_kk", randomSeed = user_seed) %>%
  visPhysics(
    stabilization = TRUE,
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(
      gravitationalConstant = -200,    # Less pull between nodes
      centralGravity = 0.01,
      springLength = 200,             # More space between connected nodes
      springConstant = 0.05,
      damping = 0.4,
      avoidOverlap = 1
    )
  ) %>%
  visLegend(
    useGroups = FALSE,
    addEdges = data.frame(
      label = c("Upregulate", "Inhibition"),
      color = c("blue", "red"),
      dashes = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ) |> 
      transform(arrows = I(list(
        list(to = list(enabled = TRUE, type = "arrow")),
        list(to = list(enabled = TRUE, type = "bar"))
      ))),
    addNodes = data.frame(
      #label = c("Max Centrality", "Active Gene", "Inactive Gene"),
      #shape = c("star", "dot", "dot"),
      #color = c("gray", "green", "gray")
      label = c("Max Centrality Per Hub"),
      shape = c("star"),
      color = c("gray")
    ),
    width = 0.2
  ) %>%
  visInteraction(
    hover = TRUE,
    dragNodes = TRUE,
    dragView = TRUE,
    zoomView = TRUE
  )
print(gg)

library(htmlwidgets)
library(webshot2)

saveWidget(gg, "network_plot.html")
webshot("network_plot.html", file = "network_plot.png", vwidth = 1600, vheight = 1200, zoom = 2)