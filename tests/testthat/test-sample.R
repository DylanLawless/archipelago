
test_that("archipelago_plot produces a ggplot and saves output files", {
  # Create dummy data
  df1 <- data.frame(set_ID = c(1,2), P = c(1, 2))
  df2 <- data.frame(P = c(0.9211610, 0.8652950, 0.6916336, 0.9191039, 0.6045032, 0.8299881),
  CHR = c( 6,  9, 12, 18, 15, 11),
  BP = c(351696, 988282, 929171, 688387, 874337, 464161),
  set_ID = c(1,2,3,4,5,6))
  
    print(df1)
    print(df2)
  
  # Use temporary file paths for output
  temp_output_path <- tempfile("test_that_file")
  temp_raw_path <-  tempfile("test_that_file")
  
  # Run the function
  p <- archipelago_plot(
    df1 = df1,
    df2 = df2,
    output_path = temp_output_path,
    output_raw = temp_raw_path
  )
  
  # Check that the returned object is a ggplot
  expect_s3_class(p, "ggplot")
  
  # Check that output files were created
  expect_true(file.exists(temp_output_path))
  expect_true(file.exists(temp_raw_path))
})
