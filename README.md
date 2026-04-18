# Modelo PK de florfenicol em tilápias do Nilo

Shiny app para simular a farmacocinética do florfenicol em *Oreochromis niloticus* a 28 °C, baseado no modelo estrutural de Tameirão et al. (2022).

> **Referência:** Tameirão ER, Rubim FM, Felix LA, Gonzaga LWF, Brandão HM,
> Murgas LDS, Ferrante M. (2022). Modelo farmacocinético de florfenicol en tilapias
> (*Oreochromis niloticus*) sometidas a diferentes temperaturas de crianza.
> *Rev Inv Vet Perú*, 33(6): e22433. DOI: [10.15381/rivep.v33i6.22433](https://doi.org/10.15381/rivep.v33i6.22433)

## Sobre o modelo

O artigo propõe um modelo de 2 compartimentos com absorção oral de primeira ordem e *lag time*. Nesta versão do app, mantive a temperatura fixa em 28 °C (os parâmetros publicados para essa condição).

Ainda não implementei o ajuste heurístico por temperatura — isso está no roadmap para uma próxima versão.

## Parâmetros do modelo (28 °C)

| Parâmetro | Valor | Unidade |
|-----------|-------|---------|
| Tlag | 7.61 | h |
| Ka | 0.83 | 1/h |
| Cl | 0.0039 | L/h/kg |
| V1 | 0.31 | L/kg |
| V2 | 0.13 | L/kg |
| Q | 0.057 | L/h/kg |

## Como usar

Instale os pacotes necessários no R:

```r
install.packages(c("shiny", "deSolve", "ggplot2", "dplyr", "bslib"))
```

Depois rode o app:

```r
shiny::runApp("app.R")
```

## Funcionalidades

- Simulação da concentração plasmática ao longo do tempo
- Comparação entre doses de 10, 15 e 20 mg/kg (valores testados no artigo)
- Tabela com os parâmetros PK publicados
- Métricas: Cmax, Cmin, AUC e tempo acima da CIM

## Sobre a CIM

O artigo não define um valor único de CIM para o gráfico. Ele cita faixas de suscetibilidade por patógeno:

- *Streptococcus agalactiae*: 0.125–16 µg/mL
- *Aeromonas* spp.: 0.125–4 µg/mL

Por isso, o app permite que o usuário defina a CIM. O valor inicial é 1 µg/mL, mas pode ser ajustado conforme o cenário.

## Limitações e próximos passos

- [ ] Ajuste por temperatura (o artigo testou 22, 28 e 32 °C)
- [ ] Intervalo de confiança para as predições (o artigo reporta CV% dos parâmetros)
- [ ] Validação com dados independentes

## Nota

Esta é uma ferramenta de apoio didático e de pesquisa. Decisões terapêuticas devem sempre ser tomadas por médico veterinário.

## Licença

MIT
