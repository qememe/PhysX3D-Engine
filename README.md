# PhysX3D Engine

Лёгкий 3D physics engine на C++20: твёрдые тела, коллизии, BVH broad-phase, импульсный солвер контактов и базовая регрессия.

## Что есть сейчас

- **Rigid body dynamics** (поступательное и вращательное движение).
- **Формы коллизий**: `Sphere`, `Box`, `Capsule`.
- **Broad-phase**: BVH (`Bounding Volume Hierarchy`).
- **Narrow-phase**:
  - sphere-sphere,
  - sphere-box,
  - box-box (SAT),
  - capsule (через аппроксимацию sphere).
- **Решение контактов**: sequential impulses + трение/реституция.
- **Сон тел** (sleeping) для снижения нагрузки.
- **Raycast** по AABB тел (приближённый broad-phase hit).
- **Проверки устойчивости**: долгие сценарии, детерминизм, NaN/Inf guard-ы.

> Проект сейчас ориентирован как foundation/ядро для дальнейшего развития (джойнты, CCD, более точные query API, улучшение параллелизма).

## Структура репозитория

```text
include/physx3d/
  math.hpp         # Vec3, Mat3, Quaternion
  collision.hpp    # Shape, AABB, Contact, collision API
  spatial.hpp      # BVH
  core.hpp         # RigidBody, PhysicsWorld
  constraints.hpp  # Base Constraint + DistanceConstraint

src/
  math.cpp
  collision.cpp
  spatial.cpp
  core.cpp
  constraints.cpp

test/
  demo.cpp         # пример симуляции (падающие кубы)
  regression.cpp   # регрессионные/стресс-тесты
```

## Требования

- CMake **3.20+**
- C++ компилятор с поддержкой **C++20** (GCC/Clang)
- Linux/macOS (проект использует `pthread`)

## Сборка

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

## Запуск

### Demo

```bash
./build/physics_demo
```

### Regression tests

```bash
ctest --test-dir build --output-on-failure
# или напрямую:
./build/physx3d_regression
```

## Быстрый пример использования

```cpp
#include <physx3d/core.hpp>
using namespace physx3d;

int main() {
    PhysicsWorld world(Vec3(0, -9.81f, 0));

    // Статика: пол
    auto* ground = world.createRigidBody();
    ground->setShape(Box(Vec3(50, 1, 50)));
    ground->setPosition(Vec3(0, -1, 0));
    ground->setStatic(true);

    // Динамика: куб
    auto* cube = world.createRigidBody();
    cube->setShape(Box(Vec3(0.5f, 0.5f, 0.5f)));
    cube->setPosition(Vec3(0, 10, 0));
    cube->setMass(1.0f);
    cube->setRestitution(0.2f);

    const float dt = 1.0f / 60.0f;
    for (int i = 0; i < 600; ++i) {
        world.step(dt);
    }

    return 0;
}
```

## Текущее покрытие регрессии

`test/regression.cpp` проверяет:

- стабильность построения BVH на наборах до 2048 объектов;
- защиту от невалидных значений (`mass`, `dt`, `NaN/Inf`);
- hit/miss кейсы raycast;
- стабильность стэка из коробок;
- long-run сценарий на сотнях тел;
- детерминизм (одинаковый input → одинаковый output);
- поведение при более крупном таймстепе.

## Ограничения (важно)

- Raycast сейчас работает по **AABB**, а не по точным surface-формам.
- `DistanceConstraint` реализован, но полноценной системы constraint management в `PhysicsWorld` пока нет.
- Для production-уровня стоит добавить:
  - continuous collision detection (CCD),
  - island solver,
  - расширенный query API,
  - сериализацию состояния мира,
  - отдельные performance benchmarks.

## Установка (опционально)

```bash
cmake --install build --prefix /usr/local
```

Это установит библиотеку в `lib` и заголовки в `include/physx3d`.

---

Если хотите, следующим шагом могу подготовить:
1) англоязычный README,
2) отдельный `CONTRIBUTING.md`,
3) минимальный benchmark suite (N тел, ms/step, memory).
