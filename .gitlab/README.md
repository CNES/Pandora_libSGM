Pour reporter les templates ou les modifier sur les autres projets Pandora. Faire:

#### preparation de la branche a pousser
```
cd le_projet_a_mettre_a_jour
git checkout master
git checkout -b template_update
git fetch origin
git rebase origin/master
```

#### 1) ajout des templates

```
git subtree add --prefix .gitlab git@gitlab.cnes.fr:OutilsCommuns/CorrelateurChaine3D/template_mr_issues.git master --squash
```

#### 2) maj des templates

```
git subtree pull --prefix .gitlab git@gitlab.cnes.fr:OutilsCommuns/CorrelateurChaine3D/template_mr_issues.git master --squash
```


