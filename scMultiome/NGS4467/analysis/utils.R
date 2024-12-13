#######################
# custom functions    #
#######################
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  #l=12500000000
  l <- format(l, scientific=TRUE)
  l <- gsub("^(\\d)e", "\\1.0e", l)
  l <- gsub("^(\\d)\\.(\\d)\\de", "\\1.\\2e", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

generate_plot_dimensions <- function(curr_msigdb_results_df, WIDTH=7.9, HEIGHT=3.4, by_variable=NULL){
    if(is.null(by_variable)){ stop("specify by_variable=<char> arg.") }
    legend_width = WIDTH * 1/8
    ID_labels_max = max(nchar(as.character(curr_msigdb_results_df$ID)))
    ID_labels_width = (ID_labels_max * (WIDTH/3) / 42)
    
    my_variable_labels_max = max(nchar(as.character(curr_msigdb_results_df[[by_variable]])))
    my_variable_width = (my_variable_labels_max * ((WIDTH / 7) / 8))
    base_width = legend_width + ID_labels_width + my_variable_width
    one_panel_x = (WIDTH / 7)
    number_of_x_panels = length(unique(curr_msigdb_results_df$comparison))
    panels_x_width = one_panel_x * number_of_x_panels
    curr_width = base_width + panels_x_width
    
    base_height = HEIGHT * 3/7
    number_of_my_variables = length(unique(curr_msigdb_results_df[[by_variable]]))
    number_of_ID_labels = 0
    for(curr_variable in unique(curr_msigdb_results_df[[by_variable]])){
        number_of_ID_labels = number_of_ID_labels + length(unique(curr_msigdb_results_df[curr_msigdb_results_df[[by_variable]] == curr_variable,]$ID))
    }
    
    panels_y_height = number_of_ID_labels * ((HEIGHT * 4/7) / 18)
    curr_height = panels_y_height + base_height
    
    return(list(curr_height, curr_width))
}

vColors = c(
    "#0000CD", "#00FF00", "#FF0000", "#808080", "#000000", "#B22222", "#DAA520",
    "#DDA0DD", "#FF00FF", "#00FFFF", "#4682B4", "#E6E6FA", "#FF8C00", "#80008B",
    "#8FBC8F", "#00BFFF", "#FFFF00", "#808000", "#FFCCCC", "#FFE5CC", "#FFFFCC",
    "#E5FFCC", "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF",
    "#FFCCFF", "#FFCCE5", "#FFFFFF", "#990000", "#666600", "#006666", "#330066",
    "#A0A0A0", "#99004C"
)

vColors20 =c(
    "#696969", "#2e8b57", "#800000", "#191970", "#808000", "#ff0000", "#ff8c00",
    "#ffd700", "#ba55d3", "#00ff7f", "#0000ff", "#adff2f", "#ff00ff", "#1e90ff",
    "#fa8072", "#dda0dd", "#ff1493", "#87cefa", "#7fffd4", "#ffe4c4"
)

vColors30 = c(
    "#000000", "#696969", "#8b4513", "#006400", "#808000", "#483d8b", "#008b8b",
    "#000080", "#9acd32", "#8fbc8f", "#800080", "#b03060", "#ff4500", "#ffa500",
    "#ffff00", "#7fff00", "#9400d3", "#00ff7f", "#dc143c", "#00ffff", "#f4a460",
    "#0000ff", "#f08080", "#da70d6", "#f0e68c", "#6495ed", "#90ee90", "#add8e6",
    "#ff1493", "#7b68ee"
)
