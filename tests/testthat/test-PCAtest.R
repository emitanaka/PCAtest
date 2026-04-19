test_that("speed test", {


  data("ex0")
  bench::mark(
    { set.seed(1)
      PCAtest(
        ex0,
        100,
        100,
        0.05,
        indload = F,
        varcorr = FALSE,
        counter = F,
        plot = F
      )
    }, {
      set.seed(1)
    PCAtest2(
      ex0,
      100,
      100,
      0.05,
      indload = F,
      varcorr = FALSE,
      counter = F,
      plot = F
    )},
    check = TRUE
  )
  
})

#   expression                                      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory     time       gc      
#  <bch:expr>                                   <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>     <list>     <list>  
#1 { set.seed(1) PCAtest(ex0, 100, 100, 0.05, … 43.6ms 43.6ms      22.9    22.1MB    229.      1    10     43.6ms <named list> <Rprofmem> <bench_tm> <tibble>
#2 { set.seed(1) PCAtest2(ex0, 100, 100, 0.05,… 34.4ms 35.4ms      28.3    19.9MB     84.9     3     9      106ms