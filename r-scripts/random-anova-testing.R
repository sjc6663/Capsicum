library(tibble)
adivA <- rownames_to_column(adivA, "ID")

res.aov <- anova_test(data = adivA1, dv = Shannon, wid = SteerID, within = Hour)

adivA %>%
  group_by(Treatment) %>%
  anova_test(Shannon ~ Hour)

adivA2 <- str_replace(adivA$SteerID, "ID", "")
adivA2

adivA1 <- adivA

adivA1$SteerID<-gsub("ID","",as.character(adivA1$SteerID))
adivA1$SteerID<-gsub("8","",as.character(adivA1$SteerID))
adivA1$SteerID<-gsub("1(\\d)","",as.character(adivA1$SteerID))
adivA1
