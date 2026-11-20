# Datos crudos de los experimentos

Esta carpeta contiene, sin procesar, los resultados de todas las corridas reportadas en el artículo, salvo la tabla de la población DTLZ2 de 1000 puntos, cuyos datos no están incluidos. El mapa entre estos archivos y las tablas del artículo está en [`../REPRODUCIBILITY.md`](../REPRODUCIBILITY.md). Las tablas se reconstruyen con `python3 scripts/make_tables.py`.

## Entorno

Según el artículo:
- Equipo: Apple M5 con 16 GB de memoria, macOS 26.6.
- Compilación: Apple clang 21 con `-O3 -march=native`.
- Límite de 10 s por cálculo.
- Instancias: generadas con pymoo 0.6.1.5 (DTLZ y WFG); las geometrías tipo MaF (lineal invertido y escalado) se generaron analíticamente.
- WFG de pygmo 2.19.8 (columna `wfg` en las carpetas 01 y 02).
- Punto de referencia: `r = z_nad + 0.1 (z_nad - z_id)`.

El driver en Python que generó las instancias y midió los tiempos **no está incluido**, y tampoco los comandos exactos con que se ejecutó.

## Carpetas

### `01_geometrias_13/`: 13 geometrías (experimento de amplitud)
- Instancias: 13 (10 frentes y 3 poblaciones), 200 puntos generados (semilla 123), d ∈ {3, 5, 8, 10, 12, 14, 16, 18, 20}. En las poblaciones, `front_size` es el número de puntos no dominados.
- Si un método supera los 10 s en una instancia, ya no se corre en dimensiones mayores de esa instancia (`skipped`).
- `fronts_comparison.csv`: `instance, kind, description, points, objectives, front_size, method, status, seconds, hv, reference_method, relative_error`.
  - Métodos: `wfg` (pygmo), `qhv`, `qhv2`, `imatrix` (bitmatrix), `wfg_imatrix` (wfg-bitmatrix).
  - `qmc` y `qmc_depth` son aproximaciones cuasi-Monte Carlo que no se reportan en el artículo.
- `fronts_wfg_c.csv`: lo mismo para WFG en C (versión 1.02 original, `wfg_c`), más `standard_error`.
- `registro_ejecucion.log`: salida en consola, en forma de tabla por instancia.
- **Alimenta:** la tabla de tiempos con d = 10 (`tab:tiempos`) y el conteo de métodos más rápidos de *Running Time by Geometry*.

### `02_tamano_maximo/`: mayor tamaño resuelto
- Frentes completos con n ∈ {10, 20, 30, 50, 75, 100, 150} y d ∈ {15, 20}.
- `tamano_maximo_frentes_llenos.txt`: tiempos de `wfg` (pygmo), `qhv2`, `imatrix` y `wfg_imatrix`. `10+` significa que superó 10 s; `---` significa que no se corrió porque ya había fallado con menos puntos.
- `tamano_maximo_wfg_original.txt`: lo mismo para `wfg_c`.
- **Alimenta:** `tab:tamano`. La columna WFG de esa tabla corresponde a `wfg_c`.

### `03_sensibilidad_componentes/`: ablación
- La misma función con componentes activados o desactivados, 200 puntos. Usa **su propia muestra de puntos**, independiente de la carpeta 01: el hipervolumen de cada instancia es distinto.
- `ablation_results.csv`: `instance, points, front_size, objectives, variant, group, flags, leaf_limit, status, seconds, nodes, leaves, hv`.
  - Variantes acumulativas: `base`, `+poda`, `+filtro`, `+orden`, `+matriz=imatrix`.
  - Variantes de bitmatrix con un componente menos: `imatrix-poda`, `imatrix-filtro`, `imatrix-orden`, `imatrix-matriz`.
  - `imatrix_ciclico`: bitmatrix con el orden cíclico de QHV-II.
- **Alimenta:** `tab:ablacion`.

### `04_estadistico_10_ejecuciones/`: estudio estadístico previo
- El mismo diseño que la carpeta 05, con 10 corridas por combinación. El artículo menciona que su mayor d resuelto coincide con el de 30 corridas.
- `statistics_results.csv`: `instance, points, objectives, seed, front_size, method, status, seconds, hv`.

### `05_estadistico_30_ejecuciones/`: estudio estadístico (resultados principales)
- Frentes cóncavo y lineal, n ∈ {10, 20, 50, 100, 200, 500, 1000}, d = 2 a 20, semillas 0 a 29 (7980 bloques).
- `statistics_results_30.csv`: las mismas columnas que la carpeta 04. Métodos: `wfg_c`, `qhv`, `qhv2`, `imatrix`, `wfg_imatrix`.
- `registro_ejecucion.log`: avance por combinación.
- **Alimenta:** `tab:maxd`, `tab:wilcoxon`, `tab:dificultad`, `tab:celdas`, `tab:mediana` y las figuras del estudio estadístico.

## Códigos de `status`

| `status` | Significado |
|---|---|
| `ok` | Terminó en menos de 10 s; `seconds` y `hv` son válidos |
| `timeout` | Superó 10 s |
| `timeout_skipped` | No se corrió: la primera semilla de la combinación superó el límite, así que las demás cuentan como fallo |
| `dominated_skip` | No se corrió: el método falló en una combinación más fácil (menos puntos y menos objetivos), así que cuenta como fallo |
| `skipped` (carpeta 01) | No se corrió: el método ya había fallado en una dimensión menor de la misma instancia |
| `n/a` o vacío (`qhv`) | No aplica: el código original de QHV admite como máximo 15 objetivos |
| `error` | El programa terminó con error |

Todo `status` distinto de `ok` cuenta como una corrida que no terminó. En las comparaciones por rangos se le asigna el último lugar del bloque (tiempo = 10 s).

## Nombres de las instancias

`concavo` (cóncavo, DTLZ2–4), `lineal` (DTLZ1), `convexo` (DTLZ2 convexo), `lineal_invertido` (MaF1/IDTLZ1), `concavo_escalado` (MaF4–5), `degenerado` (DTLZ5/MaF6), `desconectado` (DTLZ7/MaF7), `wfg1`–`wfg3`, `dtlz1_poblacion`, `dtlz2_poblacion`, `dtlz7_poblacion` (vectores de decisión aleatorios).
