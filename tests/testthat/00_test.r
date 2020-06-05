


context("Creating the constraints")


test_that("Simple order", {
  KK         <- 5
  const_list <- list( order="simple", node=1, decreasing=FALSE )
  
  const01 <- create.constraints( P1=KK, constraints=const_list )
  
  expect_is(    const01, "list" )
  expect_equal( length(const01), 6 )
  expect_equal( dim(const01$Anull), c(choose(KK,2)*2,2)  )
  
})





context("Checking function: clme")


test_that("Model fit", {
  set.seed( 42 )
  nn <- 10
  mm <- 10
  KK <- 5
  gg <- factor( rep( LETTERS[1:KK] , rep(nn,KK)), ordered=TRUE  )
  ee <- rnorm( nn*KK , 0 , 1 )
  
  xx1   <- nnet::class.ind( gg )
  xx2   <- round( runif( nn*KK , 0, 10 ), 1 )
  bb1   <- c( 20, 25, 30, 35, 40 )
  Sub   <- data.frame( Sub=rep( seq(1,mm), 5) )
  uu1a  <- data.frame( Sub=1:mm, re=rnorm( max(Sub), 0, 3 ) )
  uu1b  <- plyr::join(Sub, uu1a, by='Sub', type='left', match='all')[,-1]
  yy    <- round( xx1%*%bb1 + 1.0*xx2 + ee + uu1b , 1 )
  dat01 <- data.frame( yy, gg, xx2, Sub )
  
  const_list <- list( order="simple", node=1, decreasing=FALSE )
  out01 <- clme( yy ~ gg , data=dat01, constraints=const_list  )          
  out02 <- clme( yy ~ gg + (1|Sub), data=dat01, constraints=const_list  )  
  
  
  expect_is( out01, "clme" )
  expect_is( out02, "clme" )
  
  expect_equal( length(out01$theta), KK  )
  expect_equal( length(out02$theta), KK  )
  
  expect_equal( coef(out01), out01$theta  )
  expect_equal( coef(out02), out02$theta  )
  
  expect_is( ranef(out01), "NULL" )
  expect_is( ranef(out02), "numeric" )
  
  expect_equal( rownames( VarCorr(out01) ), c("Residual") )
  expect_equal( rownames( VarCorr(out02) ), c("Sub", "Residual") )
  
  
})












        


