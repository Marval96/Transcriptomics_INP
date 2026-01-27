# **Entornos Conda Compartidos – BMF&CL**

Este documento describe la **instalación, configuración y política de uso de Conda** en el servidor del **BMF&CL**, con el objetivo de que los entornos Conda sean **compartidos y reutilizables por todos los usuarios**, evitando entornos redundantes y promoviendo reproducibilidad.

La configuración distingue entre:

* **Entornos estables compartidos** (solo ejecución para usuarios)
* **Entornos colaborativos** (creación e instalación permitida a miembros autorizados)

---

## 1. Instalación global de Conda (root)

Conda está instalado de forma centralizada en:

```
/opt/conda
```

Para que **todos los usuarios** puedan usar Conda, se inicializa de forma global.

### 1.1 Exponer Conda a todos los usuarios

Ejecutar como **root**:

```bash
echo 'export PATH=/opt/conda/bin:$PATH' | tee /etc/profile.d/conda.sh
chmod 644 /etc/profile.d/conda.sh
```

Esto asegura que `conda` esté disponible en todas las sesiones.

### 1.2 Verificación

```bash
source /etc/profile
conda --version
```

---

## 2. Configuración recomendada de Conda

### 2.1 Evitar trabajar en el entorno `base`

Ejecutar como **root**:

```bash
conda config --system --set auto_activate_base false
```

Esto evita errores comunes y mantiene buenas prácticas.

### 2.2 Inicialización de Conda en shells de usuario

> **Nota**: si `/etc/profile.d/conda.sh` está configurado correctamente, este paso suele ser innecesario.

Cada usuario puede ejecutar **una sola vez**:

```bash
echo 'source /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc
```

---

## 3. Creación de entornos compartidos (solo ejecución)

Estos entornos son mantenidos por el administrador y **solo se ejecutan por los usuarios**.

### 3.1 Crear un entorno estable (root)

```bash
conda create -p /opt/conda/envs/lab_share python=3.10 biopython
```

Si `conda` no está disponible para root:

```bash
source /etc/profile
```

### 3.2 Activación por root (verificación)

```bash
source /opt/conda/etc/profile.d/conda.sh
conda activate /opt/conda/envs/lab_share
```

**Resultado**: entorno estable, ejecutable por todos, no editable por usuarios.

---

## 4. Entornos Conda colaborativos (editables por grupo)

Para permitir que varios usuarios **creen y mantengan entornos compartidos**, se utiliza un grupo Linux.

### 4.1 Crear grupo de trabajo

Ejecutar como **root**:

```bash
groupadd bmfclab_conda
```

### 4.2 Agregar usuarios al grupo

Ejecutar como **root**:

```bash
usermod -aG bmfclab_conda marval
usermod -aG bmfclab_conda ehatl
```

> Los usuarios deben **cerrar sesión y volver a entrar** para que el grupo tenga efecto.

Verificación (usuario):

```bash
groups
```

---

### 4.3 Preparar la ruta de entornos colaborativos

Ejecutar como **root**:

```bash
chown -R root:bmfclab_conda /opt/conda/envs
chmod 2775 /opt/conda/envs
```

* `2` (setgid): los nuevos archivos heredan el grupo
* Escritura solo permitida a miembros del grupo

---

### 4.4 Crear un entorno colaborativo (usuario miembro del grupo)

```bash
conda create -p /opt/conda/envs/lab_collab python=3.10 biopython
```

Este entorno:

* es visible para todos
* puede ser modificado por miembros de `bmfclab_conda`

### 4.5 Ajuste de permisos (solo si es necesario)

Ejecutar como **root o sudo**:

```bash
chown -R :bmfclab_conda /opt/conda/envs/lab_collab
chmod -R 2775 /opt/conda/envs/lab_collab
```

---

## 5. Buenas prácticas del laboratorio

✔ Usar entornos experimentales para pruebas

✔ Documentar cambios importantes

---

