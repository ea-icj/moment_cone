#!/usr/bin/bash

./cmp_kron_multi.py 4 6 --not_methods Orig OrigFix1
./cmp_kron_multi.py 5 4
./cmp_kron_multi.py 5 5 --not_methods Orig OrigFix1
./cmp_kron_multi.py 6 4
./cmp_kron_multi.py 6 7 --not_methods Orig OrigFix1 OrigFix2 FullProductSage TripletCache TripletCacheLoad
./cmp_kron_multi.py 7 6 --not_methods Orig OrigFix1 OrigFix2 FullProductSage
./cmp_kron_multi.py 8 3
./cmp_kron_multi.py 8 4 --not_methods Orig OrigFix1 OrigFix2 FullProductSage
./cmp_kron_multi.py 8 5 --not_methods Orig OrigFix1 OrigFix2 FullProductSage
./cmp_kron_multi.py 8 6 --timeout 30 --use_generator --methods NupletCacheMultiLevel2 NupletCacheMultiLevel3Load TripletCacheLoad
./cmp_kron_multi.py 9 5 --timeout 30 --methods NupletCache NupletCacheMultiLevel NupletCacheMultiLevel2 NupletCacheMultiLevel3Load TripletCacheLoad
./cmp_kron_multi.py 10 4 --timeout 30 --methods NupletCache NupletCacheMultiLevel NupletCacheMultiLevel2 NupletCacheMultiLevel3Load TripletCacheLoad
./cmp_kron_multi.py 12 4 --timeout 30 --use_generator --methods NupletCacheMultiLevel3Load TripletCacheLoad
