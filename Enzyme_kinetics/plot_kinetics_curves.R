library(ggplot2)
library(cowplot)

ligL_GDK <- data.frame(S = c(1, 0.5, 0.1, 0.05, 0.05, 0.025, 0.025), Vo = c(0.0022, 0.0014, 0.0005, 0.00056, 0.00057, 0.0001, 0.00009))

Vm1 = 0.0029
Km1 = 0.338

p1 <- ggplot(ligL_GDK, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm1 * x/(Km1+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 1.1) + ylim(0, 0.0025) + labs(x = "mmol/L G-diketone", y = "Reaction velocity (mmol/L*min)", title = "LigL - G-diketone")

ligL_GGE <- data.frame(S = c(0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.025, 0.025), Vo = c(0.00034, 0.0046, 0.0046, 0.0048, 0.0022, 0.0034, 0.0013, 0.0019))

Vm2 = 0.0071
Km2 = 0.055

p2 <- ggplot(ligL_GGE, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm2 * x/(Km2+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 0.15) + ylim(0, 0.006) + labs(x = "mmol/L substrate", y = "Reaction velocity (mmol/L*min)", title = "LigL - GGE")

ligL_GD <- data.frame(S = c(1.0, 0.2, 0.1, 0.025), Vo = c(0.0007, 0.0002, 0.0001, 0.00003))
Vm3 = 0.00172
Km3 = 1.53

p3 <- ggplot(ligL_GD, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm3 * x/(Km3+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 1.0) + ylim(0, 0.0007) + labs(x = "mmol/L substrate", y = "Reaction velocity (mmol/L*min)", title = "LigL - GD")

ligN_GDK <- data.frame(S = c(0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.025, 0.025), Vo = c(0.00026, 0.00025, 0.00011, 0.00022, 0.000097, 0.00095, 0.000085, 0.00008))
Vm4 = 0.0006
Km4 = 0.25

p4 <- ggplot(ligN_GDK, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm4 * x/(Km4+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 0.25) + ylim(0, 0.00027) + labs(x = "mmol/L G-diketone", y = "Reaction velocity (mmol/L*min)", title = "LigN - G-diketone")

ligN_GGE <- data.frame(S = c(0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.025, 0.025), Vo = c(0.0127, 0.014, 0.0132, 0.0107, 0.0067, 0.0074, 0.0059, 0.0042))
Vm5 = 0.018
Km5 = 0.079

p5 <- ggplot(ligN_GGE, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm5 * x/(Km5+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 0.25) + ylim(0, 0.014) + labs(x = "mmol/L G-diketone", y = "Reaction velocity (mmol/L*min)", title = "LigN - GGE")

ligD_GDK <- data.frame(S = c(0.2, 0.2, 0.1, 0.1, 0.05,  0.05, 0.025, 0.025), Vo = c(0.0003, 0.0004, 0.0002, 0.0002, 0.00009, 0.0001, 0.00003, 0.00004))
Vm6 = 0.00094
Km6 = 0.41

p6 <- ggplot(ligD_GDK, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm6 * x/(Km6+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 0.25) + ylim(0, 0.00033) + labs(x = "mmol/L G-diketone", y = "Reaction velocity (mmol/L*min)", title = "LigD - G-diketone")

ligD_GGE <- data.frame(S = c(0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.025, 0.025), Vo = c(0.0077, 0.0089, 0.0062, 0.0077, 0.0042, 0.0047, 0.0024, 0.0032))
Vm7 = 0.012
Km7 = 0.073

p7 <- ggplot(ligD_GGE, aes(x = S, y = Vo)) + geom_point() + stat_function(fun = function(x) Vm7 * x/(Km7+x)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0, 0.25) + ylim(0, 0.0095) + labs(x = "mmol/L G-diketone", y = "Reaction velocity (mmol/L*min)", title = "LigD - GGE")

plot_grid(p1, p2, p3, p4, p5, p6, p7, ncol = 3)

