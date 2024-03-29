# Pre-made R functions
source("helper_ecdf.R")

# Load pacman package for loading the other packages
if (!require("pacman")) {install.packages("pacman")}
library("pacman")

# Load required packages
p_load("data.table", "stringr", "dplyr", "tools", "parallel", "ggplot2", "terra", "rnaturalearth")

# Load all of the data
# This isn't an ideal method, but works for my purposes
load("./RData/simpsons.RData")
load("./RData/PD.RData")

ecdf_lists = get_ecdf_lists(cell_distributions)

final_ecdf_list = unlist(ecdf_lists[3:8])
all_SEs = unlist(simpsons_es[3:8])

names(final_ecdf_list) = 1:600
names(all_SEs) = 1:600
ecdf_vals = get_ecdf_vals(final_ecdf_list, all_SEs)

ntaxa = lapply(months, function(month) {
  m = as.character(month)
  return(ses.pd.results[[month]]$ntaxa)
})
ntaxa = na.exclude(unlist(ntaxa))
mean_ntaxa = mean(unlist(ntaxa))
min_ntaxa = min(unlist(ntaxa))
max_ntaxa = max(unlist(ntaxa))

names(ses.mpd.results) = months
mpd = lapply(months, function(month) {
  m = as.character(month)
  return(ses.mpd.results[[month]]$mpd.obs.z)
})
mpd = na.exclude(unlist(mpd))
mean_mpd = mean(unlist(mpd))
min_mpd = min(unlist(mpd))
max_mpd = max(unlist(mpd))

names(ses.pd.results) = months
pd = lapply(months, function(month) {
  m = as.character(month)
  return(ses.pd.results[[month]]$pd.obs.z)
})
pd = na.exclude(unlist(pd))
mean_pd = mean(unlist(pd))
min_pd = min(unlist(pd))
max_pd = max(unlist(pd))

names(ses.mntd.results) = months
mntd = lapply(months, function(month) {
  m = as.character(month)
  return(ses.mntd.results[[month]]$mntd.obs.z)
})
mntd = na.exclude(unlist(mntd))
mean_mntd = mean(unlist(mntd))
min_mntd = min(unlist(mntd))
max_mntd = max(unlist(mntd))

null_es = lapply(3:8, function(month) {
  return(simpsons_es_null[[month]])
})
null_es = unlist(null_es)

es = na.exclude(unlist(simpsons_es[3:8]))
ecdfs = na.exclude(unlist(ecdf_vals))

pct_lower_5 = length(which(ecdfs <= 0.05))/length(ecdfs)
pct_upper_5 = length(which(ecdfs >= 0.95))/length(ecdfs)
mean(ecdfs)
median(ecdfs)


# Assuming your 600 values are in a variable named data
data <- ecdfs # Example data. Replace this with your actual data.

# Create a histogram
par(mfrow=c(1,2))
hist(null_es, breaks=30, main="Null evenness values", xlab="Simpson's Evenness", ylab="Frequency", col="blue")
hist(es, breaks=30, main="Observed evenness values", xlab="Simpson's Evenness", ylab="Frequency", col="red")
qqnorm(data, main="Q-Q Plot of Data")
qqline(data, col="red")
shapiro.test(data)

mean(data)
sd(data)
min(data)
max(data)
median(data)

hist(ecdfs)

# Function to calculate densities
calculate_density <- function(data, bins) {
  # Calculate the histogram
  h <- hist(data, breaks = bins, plot = FALSE)
  
  # Convert counts to densities
  densities <- h$counts / sum(h$counts)
  
  # Return a data frame
  data.frame(mid = h$mids, density = densities, dataset = deparse(substitute(data)))
}

# Calculate densities for both datasets
bins <- seq(0.3, 1, by = 0.035)
df1_density <- calculate_density(es, bins)
df2_density <- calculate_density(null_es, bins)

# Assign dataset names
df1_density$dataset <- "Observed"
df2_density$dataset <- "Random"

# Combine into one data frame
combined_density_df <- rbind(df1_density, df2_density)

# Calculate means
mean_es <- mean(es)
mean_null_es <- mean(null_es, na.rm=TRUE)

# Create a histogram
p <- ggplot(combined_density_df, aes(x = mid, y = density, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Simpson's Evenness", y = "Percentage of All Communities", title = "Evenness by Community") +
  scale_fill_manual(name = "Dataset", values = c("Observed" = "indianred1", "Random" = "lightblue3")) +
  theme_minimal()

# Add mean lines with geom_vline
p <- p + geom_vline(xintercept = mean_es, color = "red", linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = mean_null_es, color = "blue", linetype = "dashed", size = 1.5)

# Use geom_segment to create dummy lines for the legend
p <- p + geom_segment(aes(x = Inf, y = 0, xend = Inf, yend = 0, color = "Mean Observed", linetype = "Mean Observed"), 
                      show.legend = TRUE) +
  geom_segment(aes(x = Inf, y = 0, xend = Inf, yend = 0, color = "Mean Random", linetype = "Mean Random"), 
               show.legend = TRUE)

# Adjust color and linetype scales for the dummy lines
p <- p + scale_color_manual(name = "", values = c("Mean Observed" = "red", "Mean Random" = "blue")) +
  scale_linetype_manual(name = "", values = c("Mean Observed" = "dashed", "Mean Random" = "dashed"))

# Update legends
p <- p + guides(fill = guide_legend(title = "Dataset"),
                color = guide_legend(title = ""),
                linetype = guide_legend(title = ""))

p <- p + theme(plot.title = element_text(hjust = 0.5))

# Print the plot

hist(ecdfs, breaks=50)

plot(x = na.exclude(unlist(mpd)), y = es)





# Assuming dataset1 and dataset2 are your data frames, and both have a column 'value'
# Let's add a new column to each dataset to indicate the group
es_df <- data.frame(value = es, group = "observed")
null_es_df <- data.frame(value = null_es, group = "random")

# Combine the datasets
combined_dataset <- rbind(es_df, null_es_df)

# Fit the ANOVA model
anova_model <- aov(value ~ group, data = combined_dataset)

# Display the ANOVA table
summary(anova_model)

anova_other = oneway.test(value ~ group, data = combined_dataset)

# Conducting an independent t-test
t_test_result <- t.test(value ~ group, data = combined_dataset, var.equal = FALSE)

# Print the results
print(t_test_result)

june = get_month_raster_from_tif(months=c("jun"), res=0.09)

map = ne_countries(scale='medium', type='map_units', returnclass='sf', 
                   continent=c("North America")) # Generic world map

states_provinces <- ne_states(country = c("united states of america", "canada"), returnclass = "sf")

dims = get_dims_from_tif(res=0.09)
xmin = -92.5
xmax = -67.5
ymin = 38
ymax = 48

# Example raster creation - replace with your actual raster setup
yourRaster <- rast(nrows=dims[1], ncols=dims[2]) # Adjust dimensions as necessary
ext(yourRaster) <- c(xmin, xmax, ymin, ymax) # Define your raster's extent
crs(yourRaster) <- "32618"
coordinates <- xyFromCell(yourRaster, best_indices)
points_df <- data.frame(longitude = coordinates[,1], latitude = coordinates[,2])
p <- ggplot() +
  geom_sf(data = states_provinces, fill = "darkgray", color = "lightblue") +
  geom_point(data = points_df, aes(x = longitude, y = latitude), color = "red", size = 2) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +  # Adjust xlim and ylim as needed
  theme_minimal() +
  labs(title = "Locations of 100 Study Sites") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p <- p + annotate("rect", xmin = min(points_df$longitude), xmax = max(points_df$longitude), 
                  ymin = min(points_df[10,]$latitude) - 1, ymax = max(points_df[10,]$latitude) + 1, 
                  fill = NA, color = "blue", size = 1)

bargraph.CI()

par(mfrow=c(2,2))
plot(es ~ pd, xlab="PD", ylab="Simpson's Evenness")
plot(es ~ mpd, xlab="MPD", ylab="")
plot(es ~ mntd, xlab="MNTD", ylab="Simpson's Evenness")
plot(es ~ ntaxa[-556], xlab="Species Richness", ylab="")
data = data.frame(e = es, n = ntaxa[-556])
abline(lm(data$e ~ data$n))
model <- lm(e ~ n, data = data)
summary(model)

par(mfrow=c(2,2))
plot(unlist(ecdf_vals)[-556] ~ pd, xlab="PD", ylab="ECDF(SE)")
plot(unlist(ecdf_vals)[-556] ~ mpd, xlab="MPD", ylab="")
plot(unlist(ecdf_vals)[-556] ~ mntd, xlab="MNTD", ylab="ECDF(SE)")
plot(unlist(ecdf_vals)[-556] ~ ntaxa[-556], xlab="Species Richness", ylab="")
data = data.frame(e = unlist(ecdf_vals)[-556], n = ntaxa[-556])
abline(lm(data$e ~ data$n))
model <- lm(e ~ n, data = data)
summary(model)





summed_traits_lists = mclapply(months, function(month) {
  month_pa_matrix = pa_matrices_by_month[[month]]
  trait_list = get_color_distributions(month_pa_matrix, indices=best_indices, 
                                       colors=traits_color)
  return(trait_list)
}, mc.cores=1)
names(summed_traits_lists) = months

month = "sep"
month_name = "September"

colors = c("red3", "orange3", "yellow3", "green3", "blue3", "purple3", 
           "pink3", "gray", "brown")
color_names = c("Blue", "Brown", "Gray", "Green", "Orange", "Pink", 
                "Purple", "Red", "Yellow")

row = unlist(unlist(summed_traits_lists[[month]][25,]))
names(row) = colors
row <- data.frame(Color = names(row), Count = as.numeric(row))
ggplot(row, aes(x = "", y = Count, fill = Color)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y") + 
  theme_void() +
  labs(title = month_name, fill = "Color") +
  scale_fill_manual(values = setNames(row$Color, row$Color),
                    label = color_names) +
  theme(plot.title = element_text(hjust = 0.5))


null_es_df$group_index <- (seq_along(null_es_df$value) - 1) %/% 100 + 1
mean_df <- null_es_df %>%
  group_by(group_index) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))
null_es_mean <- mean_df$mean_value[-556]

data <- data.frame(
  value = c(es, null_es_mean), # Combine the two lists
  group = factor(rep(c("Observed", "Random"), each = 599)) # Create a factor
)

bargraph.CI(x.factor = data$group, response = data$value, data = data,
            xlab = "Group", ylab = "Value", main = "Observed vs Random Values",
            col = c("blue", "red"), legend = TRUE)




# Define colors and names for the stacks
colors <- c("red3", "orange3", "yellow3", "green3", "blue3", "purple3", "pink3", "gray", "brown")
names(colors) <- c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", "Pink", "White", "Brown")

# Highly Uneven Distribution
uneven_counts <- c(7, 7, 1, 2, 2, 3, 1, 1, 1)
uneven_data <- rep(names(colors), times=uneven_counts)
uneven_df <- data.frame(Value = uneven_data, Color = uneven_data)

ggplot(uneven_df, aes(x=Value, fill=Color)) +
  geom_dotplot(binwidth=1, method="histodot", dotsize=0.5) +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.text.x=element_text(angle=45, hjust=1, size=12), 
        text=element_text(size=12),
        legend.position="none",
        plot.title=element_text(hjust=0.5),
        panel.grid.minor = element_blank()) + # Removing grid lines
  labs(fill="Stack", x="Color", y="Frequency", title="More \"Uneven\" Distribution")

# More Even Distribution
even_counts <- c(3, 2, 3, 2, 3, 3, 3, 4, 2)
even_data <- rep(names(colors), times=even_counts)
even_df <- data.frame(Value = even_data, Color = even_data)

ggplot(even_df, aes(x=Value, fill=Color)) +
  geom_dotplot(binwidth=1, method="histodot", dotsize=0.5) +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.text.x=element_text(angle=45, hjust=1, size=12), 
        text=element_text(size=12),
        legend.position="none",
        plot.title=element_text(hjust=0.5),
        panel.grid.minor = element_blank()) + # Removing grid lines
  labs(fill="Stack", x="Color", y="Frequency", title="More \"Even\" Distribution")



flowers = fread("./data/better_flowers.csv")
