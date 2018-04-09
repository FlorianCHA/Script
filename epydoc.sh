epydoc --html ./Cluster/*.py -o ./docs/ -v --graph all --inheritance grouped --docformat="restructuredtext" -N Script
git add *
git commit -m 'Mise a jour de la doc'
git push
