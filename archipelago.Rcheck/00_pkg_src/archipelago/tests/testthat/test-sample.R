test_that("archipelago_plot produces a ggplot and saves PNG output files", {
  df1 <- data.frame(set_ID = c(1,2), P = c(1, 2))
  df2 <- data.frame(P = c(0.9211610, 0.8652950, 0.6916336, 0.9191039, 0.6045032, 0.8299881),
                    CHR = c(6, 9, 12, 18, 15, 11),
                    BP = c(351696, 988282, 929171, 688387, 874337, 464161),
                    set_ID = c(1,2,3,4,5,6))
  
  temp_output_path <- tempfile("test_that_file_png")
  temp_raw_path <- tempfile("test_that_file_png")
  file_type <- "png"
  
  p <- archipelago_plot(
    df1 = df1,
    df2 = df2,
    output_path = temp_output_path,
    output_raw = temp_raw_path,
    file_type = file_type
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(paste0(temp_output_path, ".", file_type)))
  expect_true(file.exists(paste0(temp_raw_path, ".", file_type)))
})

test_that("archipelago_plot produces a ggplot and saves PDF output files", {
  df1 <- data.frame(set_ID = c(1,2), P = c(1, 2))
  df2 <- data.frame(P = c(0.9211610, 0.8652950, 0.6916336, 0.9191039, 0.6045032, 0.8299881),
                    CHR = c(6, 9, 12, 18, 15, 11),
                    BP = c(351696, 988282, 929171, 688387, 874337, 464161),
                    set_ID = c(1,2,3,4,5,6))
  
  temp_output_path <- tempfile("test_that_file_pdf")
  temp_raw_path <- tempfile("test_that_file_pdf")
  file_type <- "pdf"
  
  p <- archipelago_plot(
    df1 = df1,
    df2 = df2,
    output_path = temp_output_path,
    output_raw = temp_raw_path,
    file_type = file_type
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(paste0(temp_output_path, ".", file_type)))
  expect_true(file.exists(paste0(temp_raw_path, ".", file_type)))
})
