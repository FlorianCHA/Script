epydoc --html ./Cluster/*.py -o ./docs/ -v --graph all --inheritance grouped --docformat="restructuredtext" -n 'Script Flo'
git add *
git commit -m 'Mise a jour de la doc'
git push
