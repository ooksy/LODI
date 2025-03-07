# group variable
genogrp1 <- c("WT", "ACKG", "ACK", "HACKG", "HACK")
genogrp2 <- c("WT", "ACK/G", "HACK/G")

genogrp1 <- factor(mouse$Genotype, levels = c("WT", "ACKG", "ACK", "HACKG", "HACK"))
genogrp2 <- factor(mouse$Grouped.genotype, levels = c("WT", "ACK/G", "HACK/G"), ordered = T)

genotype1 <- sample(genogrp1, size = n, replace = T) # 우선 equal probability로 뽑음
genotype2 <- sample(genogrp2, size = n, replace = T)
