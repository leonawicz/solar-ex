library(raster)
library(parallel)
rasterOptions(chunksize=10e10, maxmemory=10e11)
path <- "/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/insolation_L48/downscaled_L48/GFDL-CM3"
files <- list.files(file.path(path, c("historical", "rcp60"), "rsds"), "tif$", full.names = TRUE)
yrs <- stringr::str_extract(files, "\\d{4}")
keep_yrs <- yrs %in% 2000:2100
files <- split(files[keep_yrs], yrs[keep_yrs])

x <- mclapply(files, function(x) calc(readAll(stack(x, quick = TRUE)), mean), mc.cores = 21) # ~220 GB peak RAM at 100 years, 32 CPUs
b <- brick(stack(x))

r <- subset(b, 1)
idx <- which(!is.na(r[]))
length(idx)

m <- extract(b, idx)

fit_slr <- function(x) as.numeric(.lm.fit(matrix(c(rep(1, 101), 2000:2100), ncol = 2), x)$coefficients)
slr <- apply(m, 1, fit_slr)

int <- slope <- r
int[idx] <- slr[1, ]
slope[idx] <- slr[2, ]
s <- stack(int, slope)

pct_change <- function(years, mod){
  int <- subset(mod, 1)
  slope <- subset(mod, 2)
  proj_fitted <- int + slope * years[2]
  clim_fitted <- int + slope * years[1]
  proj_fitted / clim_fitted - 1
}

x <- pct_change(c(2000, 2100), s)
saveRDS(x, "/workspace/UA/mfleonawicz/rsds_low48_2000-2100.rds")

# local
library(raster)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(sf)

d <- readRDS("data/rsds_ts_sample.rds") %>% mutate(id = rep(1:59, 101)) %>% group_by(id) %>% arrange(id, year) %>% mutate(crsds = cumsum(rsds))
g1 <- ggplot(filter(d, id == 1), aes(year, rsds, group = id)) + geom_line() + geom_smooth(method = "lm", color = "black", linetype = 2, se = FALSE) + labs(y = "mean annual rsds")
g2 <- ggplot(filter(d, id == 1), aes(year, crsds, group = id)) + geom_step() + labs(y = "mean annual rsds, cumulative")
ggsave("sample_cell_rsds_ts.png", g1)
ggsave("sample_cell_rsds_cts.png", g2)

#r <- 100 * readRDS("/workspace/UA/mfleonawicz/rsds_low48_2000-2100.rds")
r <- 100 * readRDS("data/rsds_low48_2000-2100.rds") # percent
#r[] <- 24 * r[] / (1000 * 0.0864) # STOP: percent change has already been computed. Do conversion upstream. # MJ/m^2/day to kWh/m^2/day
r <- aggregate(r, 100)
proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
names(r) <- "data"
pixel <- as(r, "SpatialPixelsDataFrame")
system.time( poly <- as(pixel, "SpatialPolygonsDataFrame") )
system.time( poly_sf <- as(poly, "sf") )
#us <- st_as_sf(maps::map("county", plot = FALSE, fill = TRUE))
us <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

breaks <- round(seq(min(r[], na.rm = TRUE), max(r[], na.rm = TRUE), length = 10), 2)
sfg <- scale_fill_gradient2(low = "darkblue", high = "darkred", breaks = sort(c(0, breaks[-which.min(abs(breaks))])))
gde <- guides(fill = guide_legend(title = "% \u0394", reverse = TRUE))
thm <- theme(legend.key.width = unit(2, "cm"), legend.key.height = unit(2, "cm"), legend.position = c(0.92, 0),
             panel.grid.major = element_line(colour = "transparent"))
#subtitle <- "*TO DO: Convert RSDS to another variable*" # expression("365-day mean annual Solar Irradiance"~(kWh/m^2/day))

g <- ggplot() +  geom_sf(data = poly_sf, aes(fill = data), color = "white", alpha = 0.8, size = 0.5) + 
  geom_sf(data = us, fill = "transparent", color = "#444444", size = 1) + coord_sf(crs = proj) +
  sfg + theme_map(base_size = 28) + gde + thm +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(title = "2000 - 2100 projected change in Solar Irradiance") #, subtitle = subtitle)
#ggsave("/workspace/UA/mfleonawicz/rsds_us.png", g, width = 16, height = 9.6, dpi = 300)
system.time( ggsave("rsds_us_x100.png", g, width = 2*16, height = 2*9.6, dpi = 300) )

library(tiler)
tile("rsds_us_x100.png", "data/tiles", "3-7")
