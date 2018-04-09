epydoc --html ./Cluster/*.py -o ./docs/ -v --graph all --inheritance grouped --docformat="restructuredtext" -n 'Script Flo' --navlink='<a href="https://floriancha.github.io/stage/" target="_blank">Retour au Site</a>' --top="module-tree.html"
git add *
git commit -m 'Mise a jour de la doc'
git push
