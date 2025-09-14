""" Реалізувати три алгоритми пошуку підрядка: Боєра-Мура, Кнута-Морріса-Пратта та Рабіна-Карпа.
Порівняти їхню ефективність на основі двох текстів: Стаття_1 і стаття_2"""

import timeit

# 1. Алгоритм Боєра-Мура
def boyer_moore_algo(text, pattern):
    m = len(pattern)
    n = len(text)
    if m == 0:
        return 0
    bad_char = {}
    for i in range(m):
        bad_char[pattern[i]] = i
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            return s
        else:
            s += max(1, j - bad_char.get(text[s + j], -1))
    return -1

# 2. Алгоритм Кнута-Морріса-Пратта
def kmp_algo(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    m = len(pattern)
    n = len(text)
    lps = compute_lps(pattern)
    i = j = 0
    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == m:
            return i - j
        elif i < n and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return -1

# 3. Алгоритм Рабіна-Карпа
def rabin_karp_algo(text, pattern):
    d = 256
    q = 101
    m = len(pattern)
    n = len(text)
    h = pow(d, m - 1) % q
    p = 0
    t = 0
    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q
    for s in range(n - m + 1):
        if p == t:
            if text[s:s + m] == pattern:
                return s
        if s < n - m:
            t = (d * (t - ord(text[s]) * h) + ord(text[s + m])) % q
            if t < 0:
                t += q
    return -1

with open("стаття 1.txt", encoding="utf-8") as f1:
    text1 = f1.read()
with open("стаття 2.txt", encoding="utf-8") as f2:
    text2 = f2.read()

real_substr1 = "жадібний алгоритм"
real_substr2 = "рекомендаційні системи"
fake_substr = "підрядок з іншого Всесвіту"

def measure_time(algorithm, text, pattern):
    return timeit.timeit(lambda: algorithm(text, pattern), number=10)

algorithms = {
    "Boyer-Moore": boyer_moore_algo,
    "KMP": kmp_algo,
    "Rabin-Karp": rabin_karp_algo
}
results = {}

for name, algo in algorithms.items():
    results[name] = {
        "стаття 1 (реальна))": measure_time(algo, text1, real_substr1),
        "стаття 1 (фейкова)": measure_time(algo, text1, fake_substr),
        "стаття 2 (реальна)": measure_time(algo, text2, real_substr2),
        "стаття 2 (фейкова)": measure_time(algo, text2, fake_substr),
    }

print("\nРезультати порівняння:")
for algo, timings in results.items():
    print(f"\n{algo}:")
    for case, time in timings.items():
        print(f"  {case}: {time:.6f} секунд")


