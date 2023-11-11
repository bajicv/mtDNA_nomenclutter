# Function to match cluster names/numbers across different runs
# df - dataframe with clusterings to be adjusted to each other
# cl1 - indices for the reference clustering column
# cl2 - indices for the clustering columns that will be renamed based on cl1
# cols - vector of indices of columns from which cluster names would be extracted to avoid same cluster names in different clustering runs 
# this ensures that clusters that are not present in referent column cl1 but are present in cl2 cannot have the same name as other clusters in dataset
# e.g. usage run: match_clusters(mydataframe, 4, 5, c(4,5,6,7))
# e.g. df: df <- read.csv("tmp/subset_ind_info.csv", header = T, sep = ";")

match_clusters <- function(df, cl1, cl2, cols) {
    
    max_cols <- max(df[,cols])
    
    df_tmp <- table(df[, c(cl1, cl2)]) %>%
        as.data.frame() %>%
        mutate_at(c(1,2), as.character) %>%
        mutate_at(c(1,2), as.numeric) %>%
        filter(Freq > 0)  %>%
        arrange(desc(Freq))
    
    names(df_tmp)[1:2] <- names(df)[c(cl1,cl2)]

    df_tmp$new <- NA

    # for each row in df_tmp check 
    for (r in 1:nrow(df_tmp)) {
        # if there is no value assigned in new column and if the name from ref column is not already present in new column
        if (is.na(df_tmp[r,"new"]) & !(df_tmp[r,1] %in% df_tmp$new)) {
            # get all the rows matching the current r in column and rename them according to ref column
            df_tmp[which(df_tmp[, 2] == df_tmp[r, 2]), ][, "new"] <- df_tmp[r,1]
        }
    }

    # for each row in df_tmp check 
    for (r in 1:nrow(df_tmp)) {
        # if there is no value assigned in new column and if the name from ref column is not already present in new column
        if (is.na(df_tmp[r,"new"])) {
            # get all the rows matching the current r in column and assigne them new name
            df_tmp[which(df_tmp[, 2] == df_tmp[r, 2]), ][, "new"] <- (max(max(df_tmp[, c(1,2,4)], na.rm=TRUE),max_cols) + 1)
        }
    }

    df_tmp_unique <- unique(df_tmp[,c(2,4)])  %>% mutate_if(is.factor, as.numeric)

    out_df <- suppressMessages(left_join(df, df_tmp_unique))

    out_df[, names(df)[cl2]] <- out_df$new

    out_df <- out_df  %>%
                select(-new)

    return(out_df)
}