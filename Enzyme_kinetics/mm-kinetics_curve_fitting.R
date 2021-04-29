mm <- data.frame(S = c(0.2, 0.1, 0.05, 0.025), Vo = c(0.0003, 0.0002, 0.00009, 0.00004))
model.nls <- nls(Vo ~ Vm * S/(K+S), data = mm,start = list(K = max(mm$Vo)/2, Vm = max(mm$Vo)), trace = TRUE)
summary(model.nls)

# LigD GDK - K = 0.1808706, Vm = 0.0006220
# LigD GGE - K = 0.84193, Vm = 0.04917
# LigN GDK - K = 0.08108, Vm = 0.0003354
# LigN GGE - K = 0.079431, Vm = 0.018114
# LigL GDK - K = 0.937287, Vm = 0.004211
# LigL GGE - K = 0.1475016, Vm = 0.0080681
# LigL GP - K = 1.528579, Vm = 0.001719
#
# or
# LigD GDK - K = 1.8e-4, Vm = 6.2e-7, Kcat = 0.62 (in minutes, convert to seconds)
# LigD GGE - K = 8.4e-4, Vm = 4.9e-5, Kcat = 49
# LigN GGE - K = 7.9e-5, Vm = 1.8e-5, Kcat = 18
# LigL GGE - K = 1.5e-4, Vm = 8.1e-6, Kcat = 8.1
# LigL GDK - K = 9.3e-4, Vm = 4.2e-6, Kcat = 42
# LigN GDK - K = 8.1e-5, Vm = 3.4e-7, Kcat = 3.4
# 
# Kcat = Vmax/1uM = Vmax/0.001mM = Vmax = 0.000001M
