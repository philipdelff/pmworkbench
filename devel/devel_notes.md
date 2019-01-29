# Workbench developers' notebook


## Todo
- [ ] Add a readme file for developers. Mention in .Rbuildignore

- [ ] Delete packagify branch. Philip: Done.

- [ ] Add NEWS file (for between releases)

- [x] Decide on either .Md or .html for the README file. The rest of the README files should be listed in .Rbuildignore. Philip: Done. Using .md.

- [ ] Cleanup workflows "latex_class_files", "poppk_report_latex_doNotUse". I guess we should have a structure like inst/workflows/stable inst/workflows/devel, and inst/workflows/stale or something.

I suggest shortening names. Is what is created with the templates an activity? Then I suggest adirs$scripts to be a function(...) file.path("path/to/scripts",...). Then I would use it! Ok, sure, let's discuss how to get that to work better. 

### Download_template:
- [ ] The test whether there is anything in the destination dir is unnecessarily rigid. It will be annoying in a lot of cases. I suggest creating all the destination file names and check whether any of these exist. If so, give an error and tell user which ones were found. 

- [ ] Template_available should be moved to its own file. I wouldn't expect to find other functions in download_template.R than download_template due to the name of the file. Helena: OK, will do.
Status: Helena to implement

- [ ] Definition of workflows. Guide to how to create one.
 
### List_directories
- [ ] I suggest shortening names. Is what is created with the
     templates an activity? Then I suggest adirs$scripts to be a
     function(...)  file.path("path/to/scripts",...). Then I would use
     it! Ok, sure, let's discuss how to get that to work better.
	 
- [ ] Define the list first, then do a an lapply(list,
     assign(...)). This is just to not having to check that we
     remember all of them. Make sure to define things - but only once!
     Status: philip to implement

- [ ] I don't understand what the user is supposed to use the function
for. Why not just define directories (which I suggest to call adirs or
something short in line with other names)?

### Setup03_variables.R
- [ ] Decision on whether to save setup
Philip: Should create just one object saved to an rds, not an RData. Why save all the potential crap that the user may have in their workspace? 
Helena: Yeah. I think that entire save command can actually be deleted actually. The file is short and can be sourced in the others anyway. 
- [ ] Philip: We should consider commenting out some of this. So that the user can easily include if wanted, but so that it isn't there and
     being wrong by default. If the latter is the case, we are introducing a tasks rather than help to the user. 
Helena: Ok, let's discuss 


### Todo's (to clean up?)
* Clean up s02)dataset)revew similar to done in demo


* Add testing of all functions

* All plot functions should be able to take a missingFlag and exclude that from the plots/print them as NA/missing

* generate azTheme

* Follow up on variableList for data structure function - meeting scheduled for wednesday

* What to do with the functions defined within nm templates: 
    - executeUpdate
    - cltFileUpdate
    
    
    
Add datasets:

Dataset with large range of dose groups to test the dose proportionality functions

Dataset with more than 1 study

NONMEM output file to use with r_data_staructure, nm=T

If expanding to JSON: add dataspec example for that

## Bugs

## Discussion
### General
-	Philip: I fear there is something wrong with the use of relative paths to
working dir. What is the point in not demanding a path relative to at
least QCP_MODELING? I fear that by not wanting to deal with anything
that looks like an absolute path, it is assumed that the user knows
how things should be done? What if the working dir is wrong, say from
another project that the user was working on earlier that day? I
believe we should have an internal package residing in QCP_MODELING
that gives this kind of addresses and a function that points to
addresses in our file system. Then the user must once and for all
declare a path relative to QCP_MODELING. And then, that can be changed
that one place if it should change in the future. Let's discuss, there
are some upcoming changes to the infrastructure so we will not have
the QCP_modeling structure I believe. We should have a solution that
works for the new setting. 


## Howto
- convert README.Rmd to README.md
knit("README.Rmd")