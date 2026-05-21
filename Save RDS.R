saveRDS(results_all_1, file = "results_1_new.Rds")
saveRDS(results_all_2_ext, file = "results_ext_2.Rds")
saveRDS(results_all_3_ext, file = "results_ext_3.Rds")
saveRDS(results_all_4_ext, file = "results_ext_4.Rds")
saveRDS(cph_analysis_df, file = "results_cph_CP.Rds")
saveRDS(cqr_analysis_df, file = "results_cph_CQR_width.Rds")
saveRDS(qr_analysis_df, file = "results_cph_QR_coverage.Rds")


write.table(results_all_2, file = "results_all_2.csv",
            sep = "\t", row.names = F)
write.table(results_all_3, file = "results_all_3.csv",
            sep = "\t", row.names = F)
write.table(results_all_4, file = "results_all_4.csv",
            sep = "\t", row.names = F)
write.table(results_extended_new_1, file = "results_all_1.csv",
            sep = "\t", row.names = F)


data.copy <- readRDS(file = "results_2.Rds")