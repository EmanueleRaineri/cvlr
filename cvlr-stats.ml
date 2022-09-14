open Binomlib;;

let parse_k fn =
  let f = open_in fn in
  let rec aux () =
    let l = input_line f in
    if Str.string_match (Str.regexp "^#@K:\\([0-9]+\\)") l 0 then
      let k = int_of_string (Str.matched_group 1 l)  in
      begin
        close_in f;
        k;
      end
         else aux ()
  in try aux () with End_of_file -> -1
;;

let parse_clusters f k =
  (* returns a table which associates rnames to clusters 
 and an array with the size of each cluster *)
  let tbl = Hashtbl.create 1000
  and nr = Array.init k (fun e -> 0) 
  in let rec aux () =
       let l = input_line f in
       if (l.[0] = '#')
       then
         aux ()
       else
         let fields = Array.of_list (Str.split (Str.regexp "[ \t]+") l) in
         let cl = (int_of_string fields.(1)) in
         begin
         Hashtbl.replace tbl fields.(0) cl;
         nr.(cl) <- nr.(cl) + 1;
         aux()
         end
     in try aux() with End_of_file -> (tbl , nr)
;;

let median a =
  (* 
     compute  the median of unsorted array a.
     Remove nans first  
   *)
  let nans = ref 0 in
  Array.iter (fun e -> if (Float.is_nan e ) then incr nans) a;
  (*Printf.fprintf stderr "%d nans in input array\n" !nans;*)
  let notnans = Array.length a - !nans in
  let b = Array.init notnans (fun e -> 0.0) in
  let idxb = ref 0 in
  for idxa = 0 to Array.length a - 1  do
    if (not (Float.is_nan a.(idxa))) then begin
      b.(!idxb) <- a.(idxa);
      incr idxb;
      end
  done;
  Array.sort compare b;
  let n = Array.length b   in
  match n with
  | 1 -> b.(0)
  | l when l > 1 ->
  if (0 =  l mod 2) then
    (b.( (l-1)/2 ) +. b.(l/2)) /. 2.0
  else
    b.( (l-1)/2 )
  |_ -> nan
  ;;

let gmatrix_of_file f =
  (* drnames : idx -> rname *)
  let drnames = Hashtbl.create 1000 and
      (* dgpos : idx -> gpos *)
      dgpos = Hashtbl.create 1000 and
      (* dstate: ridx, gposidx -> state *)
      dstate = Hashtbl.create 1000
      and maxridx = ref 0 and maxgposidx = ref 0 in
  let rec aux () =
    let l = input_line f in
    match ( '#' = l.[0]) with
    |true -> aux ()
    |false ->
      let fields = Array.of_list (Str.split (Str.regexp "[ \t]+") l) in
      let (ridx, rname, gposidx, gpos, state)=
        (int_of_string fields.(0),
         fields.(1),
         int_of_string fields.(2),
        int_of_string fields.(3), int_of_string fields.(4))
      in
      if (ridx > !maxridx) then begin maxridx := ridx end;
      if (gposidx > !maxgposidx) then begin maxgposidx := gposidx end;
      Hashtbl.replace drnames ridx rname;
      Hashtbl.replace dgpos gposidx gpos;
      Hashtbl.replace dstate (ridx, gposidx) state;
      ignore (aux());
      (drnames, dgpos, dstate, !maxridx, !maxgposidx)
  in try aux () with End_of_file -> (drnames, dgpos, dstate, !maxridx, !maxgposidx)
;;

let n0_n1 drnames dstate cl k d =
  (* return the number of 0s and 1s per cluster and position  *)
  let n0 = Array.init k (fun l -> Array.init d (fun j -> 0))
  and n1 = Array.init k (fun l -> Array.init d (fun j -> 0)) in
             Hashtbl.iter (fun k v -> try
                   let ridx = fst k and gposidx = snd k in
                   let clid =  Hashtbl.find cl (Hashtbl.find drnames ridx) in
                   match v with
                   | 0 -> n0.(clid).(gposidx) <- n0.(clid).(gposidx) + 1
                   | 1 -> n1.(clid).(gposidx) <- n1.(clid).(gposidx) + 1
                   | _ -> ()
                 with Not_found -> begin
                     Printf.fprintf stderr
                       "warning:read %s not found in the clustering\n"
                       (Hashtbl.find drnames (fst k));
                     ()
                   end
               ) dstate;
  (n0, n1)
;;


let abs_diff avgmeth k d =
  (* 
     compute absolute difference in methylation,
     position by position 
   *)
  let absdiff = Array.init  (k * (k-1) / 2)
                  (fun i ->  Array.init d (fun j -> 0.0)) in
  for j = 0 to d - 1 do
    let c = ref 0 in
    for l = 0 to k -1 do
      for s = l + 1 to k - 1 do
        absdiff.(!c).(j) <- Float.abs (avgmeth.(l).(j) -. avgmeth.(s).(j));
        c := !c + 1;
      done;
    done;
  done;
  absdiff
;;


let betadiff n01 n11 n02 n12 =
  Binomlib.g (n11+1) (n01+1) (n12+1) (n02+1)
;;



let methyldiff n0 n1 k d =
  let mdiff = Array.init  (k * (k-1) / 2)
                  (fun i ->  Array.init d (fun j -> 0.0)) in
  for j = 0 to d - 1 do
    let c = ref 0 in
    for l = 0 to k -1 do
      for s = l + 1 to k - 1 do
        mdiff.(!c).(j) <- betadiff n0.(l).(j) n1.(l).(j) n0.(s).(j) n1.(s).(j);
        c := !c + 1;
      done;
    done;
  done;
  mdiff
;;

let maxtbl k v1 v2 = max v1 v2
;;

let mintbl k v1 v2 = min v1 v2
;;

let _ =
  let usage_msg = "cvlr-stats <clusterfn> <matrixfn>" and
      input_files = ref [] in
    let anon_fun filename =
      input_files := !input_files@[filename]
    in let speclist = [] in
       Arg.parse speclist anon_fun  usage_msg ;
       if (List.length !input_files < 2) then
         Printf.fprintf stderr "%s\n" usage_msg
       else
         let clusterfn = List.nth !input_files 0
         and matrixfn = List.nth !input_files 1
         in 
         let k = parse_k clusterfn in
         if ( -1 = k ) then
           Printf.fprintf stderr "malformed header, missing @K\n"
         else
           let clusterf = open_in clusterfn in
           let (cl, pop) = parse_clusters clusterf k in
           let matrixf = open_in matrixfn in
           let (drnames, dgpos, dstate, maxridx, maxgposidx ) = gmatrix_of_file matrixf in
           let n = maxridx + 1  and d = maxgposidx + 1 in
           Printf.fprintf stdout "#@N:%d\n" n;
           Printf.fprintf stdout "#@D:%d\n" d;
           Printf.fprintf stdout "#@K:%d\n" k;
           Array.iteri (fun idx a -> Printf.printf "#@POP%d:%d\n" idx a) pop;
           let (n0, n1) = n0_n1  drnames dstate cl  k d in
           let avgmeth = Array.init k (fun l -> Array.init d (fun j -> 0.0))
           and totavgmeth = Array.init d (fun j -> 0.0)
           and depth = Array.init d (fun j -> 0)
           and sum0 = Array.init d (fun j -> 0)
           and sum1 = Array.init d (fun j -> 0)
           in
           for j = 0 to d -1 do
             for l = 0 to k - 1 do
               avgmeth.(l).(j) <- float_of_int n1.(l).(j) /. float_of_int (n0.(l).(j) + n1.(l).(j));
               sum0.(j) <- sum0.(j) + n0.(l).(j);
               sum1.(j) <- sum1.(j) + n1.(l).(j);
             done;
             depth.(j) <- sum0.(j) + sum1.(j);
             totavgmeth.(j) <- float_of_int sum1.(j) /. float_of_int depth.(j);
           done;
           let absdiff = abs_diff avgmeth k d and
               mdiff = methyldiff n0 n1 k d in
                      let c = ref 0 in
           for l = 0 to k -1 do
             for s = l + 1 to k - 1 do
               let adiff = absdiff.(!c) in
               let ma = median  adiff in
               Printf.printf "#@MEDIAN_ABS_DIFF_NOTNA_%d_%d:%.3f\n" l s ma ;
               incr c;
             done;
           done;
           (*  *)
           for l = 0 to k -1 do
             for s = l + 1 to k - 1 do
               Printf.printf "#@BINOMIAL_PVAL_POP%d_POP%d:%.3f\n" l s
                          (binom_pval pop.(l) pop.(s));
             done;
           done;
           let maxgpos = Hashtbl.fold maxtbl dgpos 0 and
               mingpos = Hashtbl.fold mintbl dgpos 1000000000 in
           let span = maxgpos - mingpos + 1 in
           Printf.printf "#@MINGPOS:%d\n" mingpos;
           Printf.printf "#@MAXGPOS:%d\n" maxgpos;
           Printf.printf "#@SPAN:%d\n" span;
           (* print out table *)
           for j = 0 to d - 1 do
             Printf.printf "%d\t" (Hashtbl.find dgpos j);
                 for l = 0 to k -1 do
                   Printf.printf "%d\t%d\t%.3f\t" n0.(l).(j) n1.(l).(j) avgmeth.(l).(j);
                 done;
                 Printf.printf "%d\t%d\t%.3f\t%d" sum0.(j) sum1.(j) totavgmeth.(j) depth.(j);
                 let c = ref 0 in
                 for l = 0 to k - 1 do
                   for s = l + 1 to k - 1 do
                     Printf.printf "\t%.3f" absdiff.(!c).(j);
                     Printf.printf "\t%.3g" mdiff.(!c).(j);
                     c := !c + 1;
                   done
                 done;
                 Printf.printf "\n";
           done

;;
