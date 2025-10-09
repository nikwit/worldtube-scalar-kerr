import os
import shutil
from jinja2 import Template
import sys

rel_path = str(sys.argv[1])
cwd = os.getcwd()
abs_base_path = os.path.join(cwd, rel_path)
assert os.path.exists(abs_base_path), f"Path {abs_base_path} does not exist"

ringdown_template_name = "ringdown_template.yaml"
with open(ringdown_template_name, "r") as f:
    ringdown_template = f.read()
ringdown_template = Template(ringdown_template)

ringdown_folder = os.path.join(abs_base_path, "ringdown")
os.mkdir(ringdown_folder)

ringdown_file = os.path.join(ringdown_folder, "Ringdown.yaml")

i = 0
while(os.path.exists(os.path.join(abs_base_path, f"Restart{i+1}"))):
    i+=1
volume_path = "Volume*" if i ==0 else f"Restart{i}/Volume*"

with open(ringdown_file, "w") as f:
    f.write(ringdown_template.render({"file_path": f"{str(abs_base_path)}/{volume_path}"}))

ringdown_launch_template = "/resnick/groups/sxs/nwittek/worldtube-runs/ringdown_launch_template.sh"
shutil.copy(ringdown_launch_template, os.path.join(ringdown_folder, "Submit.sh"))

os.chdir(ringdown_folder)
os.system("sbatch Submit.sh")
