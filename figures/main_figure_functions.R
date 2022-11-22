# Creating a function that plots only the relevant data from a genome-wide Manahttan plot
draw_manhattan <- function(mplots, labels, fontsize = 10, sigline_regexp = "sigline", siglabel_regexp = "siglabel") {

	stopifnot(length(labels) == length(mplots))

	# Pushing a viewport with the appropriate layout
	pushViewport(viewport(layout = grid.layout(nrow = length(mplots))))

	# Looping over all the Manhattan plots

	for(i in 1:length(mplots)) {
		# Pushing a viewport with appropriate coordinates
		pushViewport(viewport(layout.pos.row = i))
		pushViewport(plotViewport(c(0.5, 4.5, 0.5, 0.5), xscale = mplots[[i]]$vp$xscale, yscale = mplots[[i]]$vp$yscale))

		# Drawing a box around the plot
		grid.rect()

		# Plotting the points
		grid.draw(getGrob(mplots[[i]], "manhattan_points"))

		# Plotting the threshold
		grid.draw(getGrob(mplots[[i]], "manhattan_threshold"))

		# Plotting all signal lines and labels
		grid.draw(getGrob(mplots[[i]], sigline_regexp, grep = TRUE, global = TRUE))
		grid.draw(getGrob(mplots[[i]], siglabel_regexp, grep = TRUE, global = TRUE))

		# Drawing the y-axis
		grid.draw(editGrob(getGrob(mplots[[i]], "manhattan_yaxis"), gp = gpar(fontsize = fontsize)))

		# Optionally plotting the x-axis
		if(i == length(mplots)) {
			grid.draw(editGrob(getGrob(mplots[[i]], "manhattan_xaxis"), gp = gpar(fontsize = fontsize)))
			grid.draw(editGrob(getGrob(mplots[[i]], "manhattan_xlabel"), gp = gpar(fontsize = fontsize)))
		}

		# Drawing the label for that plot in the top-left corner
		grid.text(labels[i], x = 0.02, y = 0.87, just = 0)

		upViewport(2)
	}

	# Drawing the y-axis label
	upViewport()
	grid.text(expression(-log[10](italic(p))), x = grid::unit(1.2, "lines"), rot = 90, gp = gpar(fontsize = fontsize))
}

# A function that draws a set of zoomed-in Manhattan plots using common x-scale coordinates and x-axis
draw_zoomed <- function(mgrobs, labels, fontsize = 10) {

	stopifnot(length(labels) == length(mgrobs))

	vp_x <- unlist(lapply(mgrobs, function(x) x$vp$xscale))
	xscale <- range(vp_x)

	# Adjusting the font size
	for(i in 1:length(mgrobs)) {
		mgrobs[[i]] <- editGrob(mgrobs[[i]], gp = gpar(fontsize = fontsize))
	}

	pushViewport(viewport(layout = grid.layout(nrow = length(mgrobs))))

	# Plotting each grob in their own viewport
	for(i in 1:length(mgrobs)) {
		pushViewport(viewport(layout.pos.row = i))
		pushViewport(plotViewport(c(0.5, 4.5, 0.5, 0.5),
					  xscale = xscale,
					  yscale = mgrobs[[i]]$vp$yscale))
		grid.rect()
		grid.draw(getGrob(mgrobs[[i]], "pvalue_shading"))
		grid.draw(getGrob(mgrobs[[i]], "pvalue_points"))
		grid.draw(getGrob(mgrobs[[i]], "pvalue_threshold"))
		grid.draw(getGrob(mgrobs[[i]], "pvalue_feature", grep = TRUE, global = TRUE))
		grid.draw(editGrob(getGrob(mgrobs[[i]], "pvalue_yaxis"), gp = gpar(fontsize = fontsize)))

		if(i == length(mgrobs)) {
			grid.draw(editGrob(getGrob(mgrobs[[i]], "pvalue_xlabel", grep = TRUE, global = TRUE), gp = gpar(fontsize = fontsize)))
			grid.draw(editGrob(getGrob(mgrobs[[i]], "xaxis", grep = TRUE, global = TRUE), gp = gpar(fontsize = fontsize)))
		}

		# Drawing the label for that plot in the top-left corner
		grid.text(labels[i], x = 0.02, y = 0.87, just = 0)

		upViewport(2)
	}

	# Drawing the y-axis label
	upViewport()
	grid.text(expression(-log[10](italic(p))), x = grid::unit(1.2, "lines"), rot = 90, gp = gpar(fontsize = fontsize))
}

# Creating a wrapper around grid.haplotypes which crops the beginning and end of the sequence so it can be shown with a bigger font
draw_haplotypes <- function(plotting_data, difflist, plotting_range, cropping = c(0, 0), fontsize) {

	# Adjusting the plotting data based on cropping
	if(cropping[1] > 0) plotting_data <- lapply(plotting_data, function(x) x[x$pos > cropping[1], ])
	if(cropping[2] > 0) plotting_data <- lapply(plotting_data, function(x) x <- x[-((nrow(x) - cropping[2] + 1):nrow(x)), ])
	for(i in 1:length(plotting_data)) {
		plotting_data[[i]]$pos <- 1:nrow(plotting_data[[i]])
	}

	# Adjusting the difflist based on cropping
	difflist <- lapply(difflist, function(x) x <- x - cropping[1])
	difflist <- lapply(difflist, function(x) x[x %in% do.call("rbind", plotting_data)$pos])

	# Adjusting the plotting range based on cropping
	chrom <- strsplit(plotting_range, ":")[[1]][1]
	positions <- strsplit(plotting_range, ":")[[1]][2]
	start <- as.numeric(strsplit(positions, "-")[[1]][1]) + cropping[1]
	end <- as.numeric(strsplit(positions, "-")[[1]][2]) - cropping[2]

	plotting_range <- paste0(chrom, ":", as.character(start), "-", as.character(end))

	grid.haplotypes(hapdata = plotting_data, difflist = difflist, position = plotting_range, fontsize = fontsize)
}

