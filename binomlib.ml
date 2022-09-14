(*
**************************
all these routines are taken from the source code of pooled.ml  
                
***************************
*)
(*binomial table (n choose k) up to n = 20 *)
let binomial_table= 
[|
[|1.|]; 
[|1.; 1.|]; 
[|1.; 2.; 1.|]; 
[|1.; 3.; 3.; 1.|]; 
[|1.; 4.; 6.; 4.; 1.|]; 
[|1.; 5.; 10.; 10.; 5.; 1.|]; 
[|1.; 6.; 15.; 20.; 15.; 6.; 1.|]; 
[|1.; 7.; 21.; 35.; 35.; 21.; 7.; 1.|]; 
[|1.; 8.; 28.; 56.; 70.; 56.; 28.; 8.; 1.|]; 
[|1.; 9.; 36.; 84.; 126.; 126.; 84.; 36.; 9.; 1.|]; 
[|1.; 10.; 45.; 120.; 210.; 252.; 210.; 120.; 45.; 10.; 1.|];
[|1.; 11.; 55.; 165.; 330.; 462.; 462.; 330.; 165.; 55.; 11.; 1.|]; 
[|1.; 12.; 66.; 220.; 495.; 792.; 924.; 792.; 495.; 220.; 66.; 12.; 1.|]; 
[|1.; 13.; 78.; 286.; 715.; 1287.; 1716.; 1716.; 1287.; 715.; 286.; 78.; 13.; 1.|]; 
[|1.; 14.; 91.; 364.; 1001.; 2002.; 3003.; 3432.; 3003.; 2002.; 1001.; 364.; 91.; 14.; 1.|]; 
[|1.; 15.; 105.; 455.; 1365.; 3003.; 5005.; 6435.; 6435.; 5005.; 3003.; 1365.; 455.; 105.; 15.; 1.|]; 
[|1.; 16.; 120.; 560.; 1820.; 4368.; 8008.; 11440.; 12870.; 11440.; 8008.; 4368.; 1820.; 560.; 120.; 16.; 1.|]; 
[|1.; 17.; 136.; 680.; 2380.; 6188.; 12376.; 19448.; 24310.; 24310.; 19448.; 12376.; 6188.; 2380.; 680.; 136.; 17.; 1.|]; 
[|1.; 18.; 153.; 816.; 3060.; 8568.; 18564.; 31824.; 43758.; 48620.; 43758.; 31824.; 18564.; 8568.; 3060.; 816.; 153.; 18.; 1.|]; 
[|1.; 19.; 171.; 969.; 3876.; 11628.; 27132.; 50388.; 75582.; 92378.; 92378.; 75582.; 50388.; 27132.; 11628.; 3876.; 969.; 171.; 19.; 1.|]; 
[|1.; 20.; 190.; 1140.; 4845.; 15504.; 38760.; 77520.; 125970.; 167960.; 184756.; 167960.; 125970.; 77520.; 38760.; 15504.; 4845.; 1140.; 190.; 20.; 1.|]
|]
;;


let epsilon = sqrt epsilon_float
;;

let isnan (x : float) = x <> x
;;

(**
    memoizes any function of one argument
    @param f 'a->'b
    @return a'->'b
*)

let memoize f =
    let t = Hashtbl.create 1000  in 
    fun n ->
		try  Hashtbl.find t n 
    with Not_found ->    
        let res = f n  in
        (Hashtbl.add t n res;
        res) 
;;

let memoize2 f =
    let t2 = Hashtbl.create 1000
    in 
  fun n k ->
    try  Hashtbl.find t2 (n,k) 
    with Not_found ->    
        let res = f n k in
        (Hashtbl.add t2 (n,k) res;
        res) 
;;

let gammaln ( z : float ) =
    let cof =[|76.18009172947146;-86.50532032941677;
              24.01409824083091;-1.231739572450155;
              0.1208650973866179e-2;-0.5395239384953e-5|]
    in let y = ref z in
    let x = z in
    let tmp=x +.5.5 in
    let tmp = tmp -. (x +. 0.5)*.log(tmp) in
    let ser= ref 1.000000000190015 in
    for j= 0 to 5  do 
        y:=!y +. 1.0;
        ser := !ser +. cof.(j) /. !y 
    done;
-. tmp +. log(2.5066282746310005*. !ser /. x) 
;;

let factln  =
    memoize (fun (n:int) -> gammaln (float_of_int n +. 1.0)) 
;;

let logbico  = memoize2 
    (fun ( n:int ) ( k:int ) -> factln n -. factln k -. factln ( n - k) )
;;

let logpow ( e:float ) ( base:float ) = 
	if ( e = 0.0 && base = 0.0 ) then 0.0 
	else ( e *. ( log base ) )  
;;

(* binomial probability  (n choose k) p^k q^(n-k) *)

let pbico (n:int) (k:int) (p:float)  =
  let q = 1.0 -. p in
  if (n<=20) then 
    binomial_table.(n).(k) *. p ** (float_of_int k) *. q ** (float_of_int (n-k))
  else 
    exp (
        logbico n k +. 
          logpow ( float_of_int k ) p +. 
          logpow ( float_of_int n -. float_of_int k ) q
      ) 
;;


let rec sum_binom_probs n lb acc =
  if ( lb > n) then
    List.fold_left (+.) 0. acc
  else
    let p = pbico n lb 0.5 in
    sum_binom_probs n (lb+1) (p::acc)
;;

let binom_pval pop0 pop1 = 
  let s = pop0 + pop1 in
  if ( pop0 >= pop1 ) then
    sum_binom_probs s pop0 []
  else
    sum_binom_probs s pop1 []
;;


let ln_beta_function (alpha:int)  (beta:int) =
  let alpha = float_of_int alpha and beta = float_of_int beta in
  gammaln alpha +. gammaln beta -. gammaln (alpha +. beta)
;;


let h (a1:int) (b1:int) (a2:int) (b2:int) =
    let tmp = ln_beta_function (a1+a2) (b1+b2) -.
        ln_beta_function a1 b1 -.
        ln_beta_function a2 b2 in
    exp tmp
;;    

let rec g (a1:int) (b1:int) (a2:int) (b2:int) =
  match (a1,b1,a2,b2) with
  |(a1,b1,a2,b2) when (a1 = a2 && b1 = b2) -> 0.5
  |(a1,b1,a2,b2) when (a1>a2) -> g (a1 - 1) b1 a2 b2  +.
                                   (h (a1 - 1) b1 a2  b2) /. float_of_int ( a1 - 1 ) 
  |(a1,b1,a2,b2) when (a2>a1) ->  g  a1  b1 (a2 - 1) b2  -. 
                                    (h a1 b1 (a2 - 1) b2 ) /.  float_of_int ( a2 - 1 ) 
  |(a1,b1,a2,b2) when ( b1 > b2 )  ->
    g  a1  (b1 - 1)  a2 b2  -. 
      (h a1  (b1 - 1)  a2  b2) /. float_of_int  ( b1 - 1 ) 
  |(a1,b1,a2,b2) when ( b2 > b1 ) -> g  a1 b1 a2 ( b2 - 1 )  +. 
                                       ( h a1 b1 a2 ( b2 - 1 ) ) /. float_of_int   ( b2 - 1 )
  |_ -> -1.0
;;


