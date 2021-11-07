# Just insert data drame with: one column country names to plot
# and one column for values to plot


# Select map projection
# (Here it is Mollweide)


south.america <- function(title = "South America", coords = NULL) {

select_proj <- '+proj=moll'

# Select what to plot
world <- ne_countries(continent = "South America",
                      scale = "medium",
                      returnclass = "sf") %>%

         # Set map projection
         st_transform(select_proj)


# Number of countries
num <- nrow(world)


# Insert data
my_data <- data.frame(
                      name = world$name,
                      Numbers = 1 : num )


plot_data <- merge(world, my_data,
                   by = "name",
                   all.x = TRUE)

# Insert points
if (!is.null( coords ) ) {
d_points <- data.frame(long = coords[, 1] ,
                       lat  = coords[, 2]) %>%

  # Data frame to sf
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%

  # Set map projection
  st_transform(crs = select_proj)
}


##############################################
# Map theme

map_theme <- theme(

    # modify title
    plot.title = element_text(color = "chartreuse",
                              size = 16,
                              face = "bold",
                              hjust = 0.5),

    # modify axis text
    axis.title.x = element_text(color = "limegreen",
                                size = 14,
                                face = "bold.italic"),
    axis.title.y = element_text(color = "limegreen",
                                size = 14,
                                face = "bold.italic"),

    #modify background
    plot.background = element_rect(fill = "black"),
    panel.grid.major = element_line(colour = "black",
                                    size = 0.5,
                                    linetype = 3),

    panel.background = element_rect(fill = "honeydew2",
                                    colour = "honeydew2",
                                    size = 0.8,
                                    linetype = "solid"),

    # modify axis lines
    axis.line = element_line(size = 1.5,
                             colour = "white"),
    axis.ticks = element_line(size = 2,
                              colour = "white"),
    axis.text.x = element_text( size = 14,
                                colour = "white"),
    axis.text.y = element_text( size = 14,
                                colour = "white"),

    # modify legend
    legend.background = element_rect(fill = "gray",
                                     size = 0.5,
                                     linetype = "solid") )


# Admin level 0
#########################################################

# Plot

ggplot() +

  # sf object
  geom_sf(data = plot_data,
          aes(fill = my_data$Numbers)) +

   # map margins
   xlim(c(-9000000, -3470373)) +

{if (!is.null( coords ) )
  # Add points
   geom_sf(data = d_points,
           color = "black",
           size = 1,
           shape = 23)} +

  # Title
  ggtitle(title)+

  # Legend text
  labs(fill = "Number") +

  map_theme +

  scale_fill_gradientn(colors = heat.colors(num) )

#scale_fill_gradient2(low = "white",
#                     mid = "orange",
#                     high = "red",
#                     midpoint = 3)

}
