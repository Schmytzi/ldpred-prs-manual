get_scores <- function(folder){
    data <- list.files(path=folder, pattern="*profile", full.names=TRUE) %>%
        lapply(fread) %>%
        rbindlist()
    data[, sum(SCORESUM), by = IID]
}

create_model <- function(col, data){
    glm(SHBG ~ array + sex + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + col, data=data)
}

temp_predict <- function(name, models, data){
    model <- models[[name]]
    setnames(data, name, "col")
    result <- predict(model, newdata=data)
    setnames(data, "col", name)
    result
}