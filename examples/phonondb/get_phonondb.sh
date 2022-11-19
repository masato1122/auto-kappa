#
# This script downloads data of Phonondb.
# Before starting the automation claculation, users need to download Phonondb.
#
###########################
## Modify these parameters
imin=149
imax=150
###########################

for i in `seq $imin $imax`; do
    
    mpid=mp-${i}
    
    dir_name="${mpid}-20180417.tar.lzma"
    url="http://phonondb.mtl.kyoto-u.ac.jp/_downloads/${mpid}-20180417.tar.lzma"
    
    flag=`python check_url.py --url $url`
    
    if [ $flag == "False" ]; then
        continue
    fi
    
    wget $url
    
    tar xf $dir_name
    
    rm $dir_name
    
    ### You could remove this part if you want.
    dir1="${mpid}-20180417"
    dir2="${mpid}"
    mv $dir1 $dir2

done

