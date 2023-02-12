# parsevcf

## Dependencies

Requires BCFtools for an optional tag cleanup step.  
```
conda install -c bioconda bcftools
```

## Python script
```
git clone https://github.com/lilymaryam/parsevcf
cd parsevcf
python3 vcf_to_diff_script.py --help
```

## WDL script
The WDL task does not run on its own as a workflow -- it is meant to imported in a WDL workflow. It is used by the [myco](github.com/aofarrel/myco) pipeline.