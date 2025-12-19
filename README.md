# DUT Scientific Engine 3.0 - Web Visualization

Aplicação Flask para visualização interativa dos resultados do Dead Universe Theory (DUT) Engine.

## 📋 Arquivos Necessários

Certifique-se de ter os seguintes arquivos no diretório:

1. `main.py` - Aplicação Flask principal
2. `dut_engine_v3.py` - Motor científico DUT (o script original fornecido)
3. `requirements.txt` - Dependências Python

## 🚀 Deploy no Render

### Passo 1: Preparar o Repositório

1. Crie um repositório Git com os 3 arquivos acima
2. Faça commit e push para GitHub/GitLab/Bitbucket

```bash
git init
git add main.py dut_engine_v3.py requirements.txt README.md
git commit -m "Initial commit"
git push origin main
```

### Passo 2: Configurar no Render

1. Acesse [render.com](https://render.com)
2. Clique em "New +" → "Web Service"
3. Conecte seu repositório Git
4. Configure:
   - **Name**: `dut-engine` (ou outro nome)
   - **Environment**: `Python 3`
   - **Build Command**: `pip install -r requirements.txt`
   - **Start Command**: `gunicorn main:app`
   - **Instance Type**: Escolha o plano (Free tier disponível)

### Passo 3: Variáveis de Ambiente (Opcional)

Se necessário, adicione:

- `PORT`: Render define automaticamente
- `PYTHON_VERSION`: `3.11.0` (recomendado)

### Passo 4: Deploy

Clique em "Create Web Service" e aguarde o deploy (5-10 minutos).

## 💻 Execução Local

### Instalar Dependências

```bash
pip install -r requirements.txt
```

### Executar a Aplicação

```bash
python main.py
```

Acesse: `http://localhost:5000`

## 🎯 Funcionalidades

A aplicação oferece 4 modos de análise:

### 1. **Validation Mode**

- Valida parâmetros físicos do modelo DUT
- Verifica consistência de G_eff/G
- Mostra avisos de validação

### 2. **Comparison Mode** (Padrão)

- Compara DUT vs ΛCDM
- Calcula Δχ² entre os modelos
- Usa dados Pantheon+ (Supernovas)
- Aplica priors CMB (Planck 2018)

### 3. **DESI BAO Mode**

- Inclui dados de Oscilações Acústicas de Bárions (BAO)
- Usa pontos de dados DESI 2024
- Análise completa: SN + BAO + CMB

### 4. **Dimensional Test**

- Testa consistência dimensional das equações
- Verifica correção de Hdot
- Validação matemática do solver

## 📊 Visualizações

A interface exibe:

- **Cards de Resultados**: Métricas de χ² para cada modelo
- **Gráficos Interativos**: Comparação visual usando Chart.js
- **Validação Física**: G_eff/G, evolução de φ, etc.
- **Metadados**: Versão do engine, módulos habilitados, timestamp

## ⚙️ Otimizações para Produção

### 1. Cache de Resultados

O código implementa cache em memória para evitar recálculos desnecessários.

### 2. Threading

Análises pesadas rodam em threads separadas para não bloquear a interface.

### 3. Polling Assíncrono

Status updates via polling a cada 1 segundo.

### 4. Lazy Loading

Dados são carregados sob demanda (Pantheon+, covariância).

## 🔧 Troubleshooting

### Erro: "Module not found"

```bash
pip install --upgrade -r requirements.txt
```

### Timeout no Render

- Aumente o tempo limite em Settings → "Health Check Path"
- Considere usar um plano pago para mais recursos

### Memória Insuficiente

- O engine requer ~512MB RAM mínimo
- No Render Free, pode ser necessário otimizar ou usar plano pago

### Erro de Download de Dados

O engine baixa automaticamente:

- Pantheon+ SN data (~5MB)
- Covariance matrix (~50MB)

Certifique-se de que o serviço tem acesso à internet.

## 📚 Citação

Se você usar este software em pesquisa acadêmica, cite:

```
Almeida, J. (2025). Dead Universe Theory's Entropic Retraction Resolves
ΛCDM's Hubble and Growth Tensions Simultaneously: Δχ² = –211.6 with
Identical Datasets. Zenodo. https://doi.org/10.5281/zenodo.17752029
```

## 📄 Licença

Este software é disponibilizado para fins de pesquisa científica sob os termos
especificados no cabeçalho do arquivo `dut_engine_v3.py`.

## 🆘 Suporte

Para questões sobre:

- **Deploy**: Render documentation
- **DUT Engine**: ExtractoDAO Labs
- **Bugs**: Abra uma issue no repositório

---

**Versão**: 3.0
**Engine**: NASA/ESA Grade
**Status**: Production Ready ✅
