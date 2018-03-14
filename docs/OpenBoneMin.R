library(remotes)
library(dplyr)

install_github("metrumresearchgroup/OpenBoneMin")

library(OpenBoneMin)

mod <- BoneMin()


example("sim_teri")

out <- sim_teri(dose = c(20,40), dur = 9)

out %>% as_data_frame %>%count( ID)

plot(out, PTHpm ~ time|ID, scales = "same")

out <- sim_denos()

plot(out, log(DENCP) + BMDlsDENchange~.)
