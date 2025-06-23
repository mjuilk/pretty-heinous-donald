import subprocess, sys, glob
in_path = sys.argv[1]
print(in_path)
out_path = sys.argv[2]
print(out_path)
in_bams = glob.glob(f"{in_path}/*.bam")
print(len(in_bams))
counter = 0
nr = len(in_bams)

for fn in in_bams:
    fn_split = fn.split("_")[-1]
    full_out_path = f"{out_path}reads_{fn_split}"
    subprocess.run(["cp", fn, full_out_path])
    print(full_out_path)
    print(round((counter / nr) * 100, 2))
    counter += 1
