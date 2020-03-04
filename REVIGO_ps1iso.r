

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0008152","metabolic process",75.387,-0.265,-1.241, 6.986,-11.7375,0.998,0.000),
c("GO:0009653","anatomical structure morphogenesis", 1.542,-0.762, 2.793, 5.296,-5.9747,0.879,0.000),
c("GO:0009987","cellular process",63.780,-3.086, 4.800, 6.913,-13.9914,0.997,0.000),
c("GO:0033554","cellular response to stress", 2.967, 4.357,-6.211, 5.581,-7.6421,0.800,0.000),
c("GO:0048518","positive regulation of biological process", 1.744, 7.653, 1.150, 5.350,-6.9788,0.703,0.000),
c("GO:0051179","localization",18.495,-2.240, 0.566, 6.375,-6.2441,0.993,0.000),
c("GO:0065007","biological regulation",20.498,-4.698, 1.610, 6.420,-5.2700,0.993,0.000),
c("GO:0071705","nitrogen compound transport", 1.767, 1.202, 7.030, 5.355,-8.7878,0.927,0.000),
c("GO:0071840","cellular component organization or biogenesis", 8.568,-4.176, 5.591, 6.041,-7.1146,0.992,0.000),
c("GO:0006487","protein N-linked glycosylation", 0.076,-2.919,-1.132, 3.992,-5.4789,0.801,0.030),
c("GO:0008283","cell proliferation", 0.394,-3.427, 3.055, 4.704,-3.0565,0.954,0.064),
c("GO:0044237","cellular metabolic process",53.061,-3.306,-4.712, 6.833,-12.4461,0.919,0.078),
c("GO:0006807","nitrogen compound metabolic process",38.744, 3.990, 7.617, 6.696,-10.3080,0.967,0.088),
c("GO:0043170","macromolecule metabolic process",39.491,-4.715,-3.118, 6.705,-10.4498,0.929,0.089),
c("GO:0071704","organic substance metabolic process",58.357,-6.081, 1.472, 6.874,-8.9586,0.965,0.119),
c("GO:0044238","primary metabolic process",53.743,-5.529, 3.344, 6.839,-8.8601,0.965,0.120),
c("GO:0006396","RNA processing", 3.210, 0.553,-6.125, 5.615,-4.7447,0.831,0.133),
c("GO:0008219","cell death", 0.458, 5.229, 6.745, 4.769,-3.4157,0.903,0.149),
c("GO:1901360","organic cyclic compound metabolic process",30.324,-4.925,-2.511, 6.590,-3.6326,0.933,0.211),
c("GO:0055076","transition metal ion homeostasis", 0.197, 6.734, 2.657, 4.402,-4.4353,0.738,0.245),
c("GO:0044260","cellular macromolecule metabolic process",34.276,-1.767,-4.678, 6.643,-9.6696,0.832,0.249),
c("GO:0030334","regulation of cell migration", 0.144, 4.599, 3.607, 4.267,-3.8013,0.679,0.260),
c("GO:0006725","cellular aromatic compound metabolic process",29.628,-2.303,-5.772, 6.580,-4.4034,0.900,0.260),
c("GO:0046483","heterocycle metabolic process",29.664,-1.735,-6.174, 6.580,-5.2480,0.900,0.260),
c("GO:0010941","regulation of cell death", 0.344, 5.612, 0.791, 4.645,-4.3143,0.668,0.272),
c("GO:0034641","cellular nitrogen compound metabolic process",34.137,-0.550,-6.303, 6.641,-7.1073,0.868,0.277),
c("GO:0048193","Golgi vesicle transport", 0.297,-0.461, 6.709, 4.581,-4.9872,0.933,0.290),
c("GO:0015833","peptide transport", 0.298,-0.065, 6.942, 4.582,-8.1457,0.925,0.290),
c("GO:0018196","peptidyl-asparagine modification", 0.015,-5.329,-0.625, 3.291,-3.2373,0.871,0.296),
c("GO:0010646","regulation of cell communication", 0.929, 7.467,-0.230, 5.076,-3.8861,0.677,0.303),
c("GO:0023051","regulation of signaling", 0.934, 7.474,-0.756, 5.079,-3.9830,0.691,0.303),
c("GO:0010604","positive regulation of macromolecule metabolic process", 0.988, 5.842,-0.834, 5.103,-6.7645,0.556,0.305),
c("GO:0050790","regulation of catalytic activity", 1.575, 7.321, 1.395, 5.306,-4.9066,0.656,0.307),
c("GO:0065009","regulation of molecular function", 1.726, 7.319, 0.702, 5.345,-5.2434,0.715,0.311),
c("GO:1905244","regulation of modification of synaptic structure", 0.000, 3.653, 2.015, 0.903,-3.1746,0.796,0.316),
c("GO:0006109","regulation of carbohydrate metabolic process", 0.090, 6.229, 1.632, 4.061,-3.2426,0.716,0.322),
c("GO:0048523","negative regulation of cellular process", 1.830, 7.112, 1.654, 5.371,-5.9318,0.642,0.329),
c("GO:0043412","macromolecule modification", 9.785,-0.297,-4.269, 6.099,-4.1931,0.893,0.331),
c("GO:0051128","regulation of cellular component organization", 1.586, 6.396, 0.544, 5.308,-5.5735,0.647,0.333),
c("GO:0048519","negative regulation of biological process", 1.984, 7.562, 0.387, 5.406,-5.1475,0.700,0.334),
c("GO:0016192","vesicle-mediated transport", 1.085, 0.475, 7.196, 5.144,-5.7932,0.930,0.334),
c("GO:0009894","regulation of catabolic process", 0.146, 6.826,-0.824, 4.272,-5.5058,0.696,0.336),
c("GO:0051641","cellular localization", 2.041, 0.551, 6.572, 5.418,-5.0841,0.927,0.347),
c("GO:0046907","intracellular transport", 1.564, 0.888, 7.102, 5.302,-7.3279,0.902,0.349),
c("GO:1990928","response to amino acid starvation", 0.000, 3.306,-5.812, 1.602,-3.6990,0.897,0.363),
c("GO:1901796","regulation of signal transduction by p53 class mediator", 0.012, 6.019,-2.949, 3.189,-3.5100,0.674,0.367),
c("GO:1901698","response to nitrogen compound", 0.178, 4.559,-6.533, 4.359,-4.0926,0.880,0.375),
c("GO:0034248","regulation of cellular amide metabolic process", 0.700, 6.268,-1.040, 4.954,-6.2277,0.647,0.391),
c("GO:0098586","cellular response to virus", 0.013, 3.743,-6.175, 3.236,-3.5100,0.904,0.391));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
head(one.data);

one.data[order(one.data$log10_p_value), ]


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
